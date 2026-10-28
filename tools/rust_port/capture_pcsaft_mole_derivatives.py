"""Capture PC-SAFT mole gradients from the prebuilt exact-source oracle.

Each case runs in an isolated, time-bounded process. The oracle's F_n values
are saved directly, without using native model outputs as expected values.
The helper evaluates the original nonassociation kernels directly; it does not
include association contributions from the public model dispatch. Build the
helper against the same source DLL and module files, in an isolated directory:

    gfortran -O3 -fdefault-real-8 -fdefault-double-8 -frecursive -shared \
      -I <source>/.reference tools/rust_port/pcsaft_mole_probe.f90 \
      <source>/.reference/thermopack.dll -o .reference/pcsaft-mole/pcsaft_mole_probe.dll

The source DLL and its modules are read-only; use the compiler matching its
metadata. Failed source calls are recorded in the fixture metadata.
"""
import argparse, csv, hashlib, io, json, math, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = [
    ("CO2,C1", 300., 1e-3, [0.05, 0.05]),
    ("CO2,C1", 320., 5e-4, [0.1, 0.2]),
    ("N2,CO2", 280., 2e-4, [0.1, 0.1]),
    ("C1,C2,C3", 350., 1e-3, [0.02, 0.03, 0.05]),
    ("CO2,C1,C3", 320., 5e-5, [0., 0.04, 0.06]),
    ("CO2,C1,C3", 400., 2e-4, [0.097, 0.002, 0.001]),
    ("AR,KR,N2,O2", 330., 1e-3, [0.01, 0.02, 0.03, 0.04]),
    ("C2_1,C2,C3,IC4", 400., 8e-4, [0., 0.02, 0.03, 0.05]),
    ("NC4,NC5,NC6,NC7", 450., 1e-3, [0.01, 0.02, 0.03, 0.04]),
    ("NC8,NC9,NC10,NC11", 550., 1e-3, [0.01, 0.02, 0.03, 0.04]),
    ("NC12,NC13,NC14,NC15", 600., 1e-3, [0.01, 0.02, 0.03, 0.04]),
    ("NC16,NC17,NC18,NC19", 650., 1e-3, [0.01, 0.02, 0.03, 0.04]),
    ("NC20,NC21,NC22,NC23", 700., 1e-3, [0.01, 0.02, 0.03, 0.04]),
    ("NC24,NC25,IC5", 700., 1e-3, [0.02, 0.03, 0.05]),
    ("BENZENE,TOLU,CYCLOHEX", 450., 1e-3, [0.02, 0.03, 0.05]),
    ("CO2,H2S,DME", 350., 1e-3, [0.02, 0.03, 0.05]),
    ("C2_1,ACETYLENE,CO2", 350., 1e-3, [0.02, 0.03, 0.05], "Gross2005"),
    ("C2_1,ACETYLENE,BENZENE", 400., 1e-3, [0.0, 0.04, 0.06], "Gross2005ADJQ"),
    ("CO2,C1", 350., 5e-4, [0.04, 0.06], "Tang_Gross2010"),
]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record_cases():
    text = (ROOT / "rust/src/pc_saft_data.rs").read_text(encoding="utf-8")
    pure = text.split("pub(super) const PURE_RECORDS:", 1)[1].split("pub(super) const", 1)[0]
    pattern = (r'PcSaftRecord\s*\{\s*ident: "([^"]+)"\s*,\s*m: ([^,]+),'
               r'\s*sigma: ([^,]+),\s*eps_depth_divk: ([^,]+),\s*reference: "([^"]+)"')
    records = re.findall(pattern, pure, re.S)
    cases = []
    for ident, m, sigma, eps, reference in records:
        other = "CO2" if ident == "C1" else "C1"
        # The second component uses its default fit. The target fit determines
        # its own diameter and the chosen volume gives eta=0.325 for the pair.
        other_m, other_sigma, other_eps = ((2.0729, 2.7852e-10, 169.21)
            if other == "CO2" else (1.0, 3.7039e-10, 150.03))
        t = 380.0
        d = float(sigma) * (1.0 - 0.12 * math.exp(-3.0 * float(eps) / t))
        od = other_sigma * (1.0 - 0.12 * math.exp(-3.0 * other_eps / t))
        dense_volume = math.pi * 6.02214076e23 / (6.0 * 0.325) * (
            0.04 * float(m) * d**3 + 0.06 * other_m * od**3)
        cases.append((f"{ident},{other}", t, 4e-4, [0.04, 0.06], reference, "record_dilute"))
        cases.append((f"{ident},{other}", t, dense_volume, [0.04, 0.06], reference, "record_dense"))
    if len(records) != 77:
        raise RuntimeError(f"expected 77 records, found {len(records)}")
    return cases


def audit_data(metadata, archive, fixture_text, historical, source_hashes):
    """Validate the complete capture in memory, independently of stored labels."""
    def require(condition, message):
        if not condition:
            raise RuntimeError(message)

    hash_fields = ("fixture_sha256", "archive_sha256", "library_sha256",
                   "source_metadata_sha256", "verified_library_sha256", "helper_compiler_sha256",
                   "probe_sha256", "probe_source_sha256", "worker_source_sha256")
    for key in hash_fields:
        require(isinstance(metadata.get(key), str) and re.fullmatch(r"[0-9a-f]{64}", metadata[key]),
                f"invalid provenance hash: {key}")
    for group in ("sources", "helper_modules_sha256"):
        entries = metadata.get(group)
        require(isinstance(entries, dict) and bool(entries), f"invalid provenance: {group}")
        for value in entries.values():
            require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value),
                    f"invalid provenance hash in {group}")
    require(set(metadata["helper_modules_sha256"]) == {"thermopack_var.mod", "pc_saft_nonassoc.mod"},
            "unexpected helper modules")
    for key in ("baseline", "compiler", "flags", "lapack_commit", "library_sha256", "sources"):
        require(metadata.get(key) == historical.get(key) and key in historical,
                f"historical metadata mismatch: {key}")
    require(metadata["verified_library_sha256"] == metadata["library_sha256"], "source DLL mismatch")
    require(metadata.get("helper_compiler") == metadata["compiler"], "helper compiler mismatch")
    require(metadata.get("helper_flags") == ["-O3", "-fdefault-real-8", "-fdefault-double-8", "-frecursive", "-shared"], "helper flags mismatch")
    for key, value in source_hashes.items():
        require(metadata[key] == value, f"source provenance mismatch: {key}")
    require(metadata.get("stored_records") == 77, "record count mismatch")

    expected = []
    for case in CASES + record_cases():
        components, temperature, volume, moles = case[:4]
        for simplified in (False, True):
            expected.append(dict(mode="PCSAFT_STATE", components=components, temperature=temperature,
                volume=volume, moles=moles, simplified=simplified,
                reference=case[4] if len(case) > 4 else "Default",
                category=case[5] if len(case) > 5 else "targeted"))
    require(len(archive) == len(expected) == 346, "archive count mismatch")
    require([item.get("request") for item in archive] == expected, "archive requests do not cover the declared case matrix")
    fixture = io.StringIO(newline="")
    writer = csv.writer(fixture, delimiter=";", lineterminator="\n")
    writer.writerow(["components", "temperature", "volume", "moles", "f_n", "simplified", "reference"])
    successes = failures = 0
    for item in archive:
        rc, stdout, stderr = item.get("returncode"), item.get("stdout"), item.get("stderr")
        require(type(rc) is int and isinstance(stdout, str) and isinstance(stderr, str),
                "censored or malformed source process")
        lines = [line for line in stdout.splitlines() if line.startswith("RESULT=")]
        if rc != 0:
            require(not lines, "failed process emitted a result")
            classification = "failure"
            failures += 1
        else:
            require(len(lines) == 1, "successful process lacks a unique raw result")
            try:
                values = json.loads(lines[0].split("=", 1)[1])
            except (ValueError, TypeError) as error:
                raise RuntimeError("malformed raw source result") from error
            request = item["request"]
            require(isinstance(values, list) and len(values) == len(request["moles"]) and
                    all(type(x) in (int, float) and math.isfinite(x) for x in values),
                    "invalid raw source result")
            classification = "success"
            successes += 1
            writer.writerow([request["components"], repr(request["temperature"]), repr(request["volume"]),
                ",".join(map(repr, request["moles"])), ",".join(format(x, ".17e") for x in values),
                str(request["simplified"]).lower(), request["reference"]])
        require(item.get("outcome") == classification, "stored outcome disagrees with raw process")
    require(fixture.getvalue() == fixture_text, "fixture does not reconstruct from raw successful source output")
    require(metadata.get("attempted_cases") == len(archive), "attempted count mismatch")
    require(metadata.get("cases") == successes, "successful count mismatch")
    require(metadata.get("failure_count") == failures, "failure count mismatch")
    return len(archive), successes, failures


def audit():
    fixture_dir = ROOT / "tests/fixtures"
    meta_path = fixture_dir / "pcsaft_mole_derivatives.metadata.json"
    metadata = json.loads(meta_path.read_text(encoding="utf-8"))
    for key, name in [("fixture_sha256", "pcsaft_mole_derivatives.csv"),
                      ("archive_sha256", "pcsaft_mole_derivatives.oracle.jsonl")]:
        actual = digest(fixture_dir / name)
        if actual != metadata.get(key):
            raise RuntimeError(f"{name}: expected {metadata.get(key)}, got {actual}")
    archive = [json.loads(line) for line in (fixture_dir / "pcsaft_mole_derivatives.oracle.jsonl").read_text(encoding="utf-8").splitlines()]
    historical = json.loads((fixture_dir / "pcsaft_state.metadata.json").read_text(encoding="utf-8"))
    source_hashes = {key: digest(ROOT / "tools/rust_port" / name) for key, name in [
        ("probe_source_sha256", "pcsaft_mole_probe.f90"),
        ("worker_source_sha256", "pcsaft_mole_oracle.py")]}
    counts = audit_data(metadata, archive, (fixture_dir / "pcsaft_mole_derivatives.csv").read_text(encoding="utf-8"), historical, source_hashes)
    print(f"offline audit: {counts[0]} calls, {counts[1]} successes, {counts[2]} retained failures; fixture reconstructed from raw results")


def main(source, compiler_directory):
    source = source.resolve()
    compiler_directory = compiler_directory.resolve()
    reference_dir = source / ".reference"
    metadata = json.loads((reference_dir / "metadata.json").read_text(encoding="utf-8"))
    library = reference_dir / "thermopack.dll"
    if digest(library) != metadata["library_sha256"]:
        raise RuntimeError("source DLL hash does not match source metadata")
    compiler = compiler_directory / "gfortran.exe"
    version = subprocess.check_output([str(compiler), "--version"], text=True, encoding="utf-8").splitlines()[0]
    if version != metadata["compiler"]:
        raise RuntimeError(f"compiler does not match source metadata: {version}")
    output = ROOT / ".reference/pcsaft-mole"
    output.mkdir(parents=True, exist_ok=True)
    helper_source = ROOT / "tools/rust_port/pcsaft_mole_probe.f90"
    helper = output / "pcsaft_mole_probe.dll"
    flags = ["-O3", "-fdefault-real-8", "-fdefault-double-8", "-frecursive", "-shared"]
    build = [str(compiler), *flags, "-I", str(reference_dir), str(helper_source), str(library), "-o", str(helper)]
    subprocess.run(build, check=True, capture_output=True, text=True, encoding="utf-8", timeout=120)
    worker = ROOT / "tools/rust_port/pcsaft_mole_oracle.py"
    results, outcomes = [], []
    all_cases = CASES + record_cases()
    for case in all_cases:
        components, temperature, volume, moles = case[:4]
        reference = case[4] if len(case) > 4 else "Default"
        category = case[5] if len(case) > 5 else "targeted"
        for simplified in (False, True):
            request = dict(mode="PCSAFT_STATE", components=components, temperature=temperature,
                           volume=volume, moles=moles, simplified=simplified, reference=reference,
                           category=category)
            command = [sys.executable, str(worker), str(source), str(compiler_directory)]
            try:
                process = subprocess.run(command, input=json.dumps(request)+"\n", text=True,
                    encoding="utf-8", errors="replace", capture_output=True, timeout=45)
                item = dict(request=request, returncode=process.returncode, stdout=process.stdout,
                            stderr=process.stderr, outcome="failure")
                lines = [line for line in process.stdout.splitlines() if line.startswith("RESULT=")]
                if process.returncode or len(lines) != 1:
                    item["reason"] = "source process failed or omitted a unique result"
                else:
                    response = json.loads(lines[0].split("=", 1)[1])
                    if len(response) != len(moles) or not all(math.isfinite(x) for x in response):
                        item["reason"] = "wrong result count or nonfinite source value"
                    else:
                        item["outcome"] = "success"
                        results.append((components, temperature, volume, moles, response, simplified, reference))
            except subprocess.TimeoutExpired as error:
                def decode(value):
                    return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value or ""
                item = dict(request=request, returncode=None, stdout=decode(error.stdout),
                            stderr=decode(error.stderr), outcome="timeout", reason="45 second limit")
            except (ValueError, TypeError) as error:
                item["reason"] = f"invalid source output: {error}"
            outcomes.append(item)
    fixture_dir = ROOT / "tests/fixtures"
    path = fixture_dir / "pcsaft_mole_derivatives.csv"
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter=";", lineterminator="\n")
        writer.writerow(["components", "temperature", "volume", "moles", "f_n", "simplified", "reference"])
        for components, t, v, n, values, simplified, reference in results:
            writer.writerow([components, repr(t), repr(v), ",".join(map(repr, n)),
                             ",".join(format(x, ".17e") for x in values), str(simplified).lower(), reference])
    archive = fixture_dir / "pcsaft_mole_derivatives.oracle.jsonl"
    archive.write_text("".join(json.dumps(item, ensure_ascii=True, sort_keys=True, separators=(",", ":")) + "\n"
                               for item in outcomes), encoding="utf-8", newline="\n")
    modules = ["thermopack_var.mod", "pc_saft_nonassoc.mod"]
    metadata.update(generator="tools/rust_port/capture_pcsaft_mole_derivatives.py",
        fixture_sha256=digest(path), archive_sha256=digest(archive), cases=len(results),
        attempted_cases=len(outcomes), stored_records=len(record_cases()) // 2,
        failure_count=sum(item["outcome"] != "success" for item in outcomes),
        source_metadata_sha256=digest(reference_dir / "metadata.json"),
        verified_library_sha256=digest(library), helper_compiler=version,
        helper_compiler_sha256=digest(compiler), helper_flags=flags,
        helper_modules_sha256={name: digest(reference_dir / name) for name in modules},
        probe_sha256=digest(helper), probe_source_sha256=digest(helper_source),
        worker_source_sha256=digest(worker), convention="direct nonassociation F_PC_SAFT_TVn or F_sPC_SAFT_TVn F_n at fixed temperature and volume")
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(f"captured {len(results)}/{len(outcomes)} calls; all outcomes retained")
    audit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path)
    parser.add_argument("--compiler-directory", type=Path)
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.audit:
        audit()
    elif args.source and args.compiler_directory:
        main(args.source, args.compiler_directory)
    else:
        parser.error("provide --source and --compiler-directory, or --audit")
