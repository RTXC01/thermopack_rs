"""Compare Rust with a fresh executable linked to the exact-source Fortran build."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/fixtures/saftvrmie_tvn.csv"
PROBE = Path(__file__).with_name("saftvrmie_probe.f90")


def cases():
    def case(comps, t, v, amounts, order=-1, hs=1, terms=31, cutoff=0, shift=0, kind=0, ref="DEFAULT", variant="default"):
        return [comps, ref, order, hs, terms, cutoff, shift, kind, t, v, amounts, variant]
    # Extensive TVN residual Helmholtz (kind=0)
    yield case("CO2,C1", 280, 0.0004, "0.8,1.2")
    yield case("NH3", 320, 0.0003, "2")
    yield case("H2O,MEOH", 380, 0.0003, "0.8,1.2")
    yield case("ETOH,C1", 350, 0.0004, "1.2,0.8")
    for hs in range(2, 6):
        yield case("CO2,C1", 280, 0.0004, "0.8,1.2", hs=hs)
    for terms in (1, 2, 4, 8, 16):
        yield case("CO2,C1", 280, 0.0004, "0.8,1.2", terms=terms)
    for shift in (0, 1):
        yield case("CO2", 300, 0.0004, "2", cutoff=3.5, shift=shift)
    for order in (0, 1, 2):
        yield case("H2", 40, 0.0005, "2", order=order, hs=4)
        yield case("H2,NE", 54, 0.0005, "1.3,0.7", order=order, hs=5)
    # DFT molar dispersion (kind=1)
    yield case("H2O,MEOH", 380, 0.0003, "0.8,1.2", kind=1)
    yield case("H2,HE", 65, 0.0007, "0.7,1.3", order=2, hs=4, kind=1)
    yield case("CO2", 300, 0.0004, "2", cutoff=3.5, shift=1, kind=1)
    yield case("H2", 45, 0.0005, "2", order=2, hs=4, cutoff=3.5, shift=1, kind=1)
    # Dufal association RDF (kind=2)
    yield case("H2O,MEOH", 380, 0.0003, "0.8,1.2", kind=2)
    # Raw non-unit mole scaling & Elliot combining rule (kind=0)
    yield case("H2O,MEOH", 380, 0.0015, "4,6")
    yield case("H2O,MEOH", 380, 0.0003, "0.8,1.2", variant="elliot")
    yield case("H2,NE", 54, 0.0005, "1.3,0.7", order=1, hs=5, variant="additive-exact-quantum")
    # Independent submodel options (kind=0)
    yield case("H2", 45, 0.0005, "2", order=2, hs=4, cutoff=3.5, shift=1, variant="jaen-kahn")
    for variant in ("geometric", "arithmetic-cross", "linear", "bmcsl", "kolafa",
                    "pair-a3", "classical-a2"):
        yield case("H2,NE", 54, 0.0005, "1.3,0.7", order=2, hs=4, variant=variant)
    yield case("H2,NE", 54, 0.00025, "0.65,0.35", order=2, hs=4, variant="linear-then-nonadd")
    yield case("CO2,C1", 280, 0.0004, "0.8,1.2", hs=5, variant="no-extra")
    # DFT molar hard-sphere (kind=3)
    yield case("H2O,MEOH", 380, 0.0003, "0.8,1.2", kind=3)
    yield case("H2,HE", 65, 0.0007, "0.7,1.3", order=2, hs=4, kind=3)
    yield case("CO2,C1", 280, 0.0004, "0.8,1.2", hs=5, kind=3, variant="no-extra")
    yield case("CO2,C1", 280, 0.0004, "0.8,1.2", hs=5, kind=3)
    # Normal TP path: v = P [Pa], shift = phase (1=LIQ, 2=VAP) (kind=4)
    yield case("CO2,C1", 280, 5.0e6, "0.4,0.6", shift=2, kind=4)
    yield case("CO2,C1", 250, 8.0e6, "0.8,0.2", shift=1, kind=4)
    yield case("NH3", 320, 5.0e6, "1", shift=1, kind=4)
    yield case("H2O,MEOH", 380, 5.0e6, "0.4,0.6", shift=1, kind=4)
    yield case("H2O,MEOH", 380, 5.0e6, "0.4,0.6", shift=1, kind=4, variant="elliot")
    yield case("ETOH,C1", 350, 3.0e6, "0.6,0.4", shift=2, kind=4)
    yield case("H2", 40, 2.0e6, "1", order=1, hs=4, shift=2, kind=4)
    yield case("H2,NE", 54, 3.0e6, "0.65,0.35", order=2, hs=4, shift=2, kind=4)
    # Model-side diameters, effective parameters, alpha, de Broglie (kind=5)
    yield case("CO2,C1", 280, 0.0004, "0.4,0.6", kind=5)
    yield case("CO2", 300, 0.0004, "1", cutoff=3.5, shift=1, kind=5)
    yield case("CO2", 300, 0.0004, "1", cutoff=1000, shift=1, kind=5)
    yield case("H2,NE", 54, 0.0005, "0.65,0.35", order=2, hs=4, kind=5)
    yield case("H2,NE", 54, 0.0005, "0.65,0.35", order=2, hs=4, kind=5, variant="arithmetic-cross")
    # Pure-fluid saturation pressure and coexisting roots (kind=6)
    yield case("H2O", 300, 0, "1", kind=6)
    yield case("H2O", 325, 0, "1", kind=6)
    yield case("H2O", 373.15, 0, "1", kind=6)
    yield case("H2O", 679.0, 0, "1", kind=6)
    yield case("H2O", 679.03, 0, "1", kind=6)
    yield case("CO2", 260, 0, "1", kind=6)
    # Pair potential (default 500*eps cap) and uncapped force: v = r / sigma_ij (kind=7)
    yield case("CO2", 260, 0.5, "1", kind=7)
    yield case("CO2", 260, 1.3, "1", kind=7)
    yield case("CO2", 300, 1.3, "1", cutoff=3.5, shift=1, kind=7)
    yield case("CO2", 300, 3.6, "1", cutoff=3.5, shift=1, kind=7)
    yield case("H2,NE", 54, 1.2, "0.65,0.35", order=2, hs=4, kind=7)


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--write-fixture", action="store_true")
    args = parser.parse_args()
    ref = ROOT / ".reference"
    compiler = Path(shutil.which(args.compiler) or args.compiler).resolve()
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([str(ref), str(compiler.parent), env["PATH"]])
    lib = ref / ("thermopack.dll" if os.name == "nt" else "libthermopack.so")
    flags = [] if os.name == "nt" else ["-Wl,-rpath," + str(ref)]
    all_cases = list(cases())
    with tempfile.TemporaryDirectory(dir=ref) as tmp:
        exe = Path(tmp) / ("mie.exe" if os.name == "nt" else "mie")
        subprocess.run([compiler, "-fdefault-real-8", "-fdefault-double-8", "-I", ref,
                        PROBE, lib, *flags, "-o", exe], check=True, env=env)
        live = Path(tmp) / "live.csv"
        with live.open("w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file, delimiter=";", lineterminator="\n")
            writer.writerow("components ref order hs terms cutoff shift dispersion T V n variant blocks dv_fv dt_fv".split())
            for c in all_cases:
                comps, param, order, hs, terms, cutoff, shift, kind, t, v, amounts, variant = c
                nc = len(comps.split(","))
                stdin = f"{comps}\n{param}\n{variant}\n{nc} {order} {hs} {terms} {cutoff} {shift} {kind}\n{t} {v}\n{amounts}\n"
                proc = subprocess.run([exe], input=stdin, text=True, env=env, capture_output=True,
                                      timeout=30)
                if proc.returncode:
                    raise RuntimeError((c, proc.stdout, proc.stderr))
                if kind == 4:
                    blocks, = [s.split()[1:] for s in proc.stdout.splitlines() if s.startswith("TP ")]
                    expected_len = 16 + 8 * nc + nc * nc
                    if len(blocks) != expected_len or not all(math.isfinite(float(val)) for val in blocks):
                        raise ValueError((c, len(blocks), expected_len, blocks))
                    writer.writerow([*c, ",".join(blocks), "0", "0"])
                elif kind == 5:
                    blocks, = [s.split()[1:] for s in proc.stdout.splitlines() if s.startswith("HOOKS ")]
                    expected_len = 5 + 3 * nc + nc * nc
                    if len(blocks) != expected_len:
                        raise ValueError((c, len(blocks), expected_len, blocks))
                    numeric = blocks[:-1] if order < 0 else blocks
                    if not all(math.isfinite(float(val)) for val in numeric):
                        raise ValueError((c, blocks))
                    if order < 0 and blocks[-1] != "UNSUPPORTED_SOURCE":
                        raise ValueError((c, "expected UNSUPPORTED_SOURCE marker", blocks[-1]))
                    writer.writerow([*c, ",".join(blocks), "0", "0"])
                elif kind == 6:
                    blocks, = [s.split()[1:] for s in proc.stdout.splitlines() if s.startswith("SAT ")]
                    if len(blocks) != 3 or not all(math.isfinite(float(val)) and float(val) > 0.0 for val in blocks):
                        raise ValueError((c, blocks))
                    writer.writerow([*c, ",".join(blocks), "0", "0"])
                elif kind == 7:
                    blocks, = [s.split()[1:] for s in proc.stdout.splitlines() if s.startswith("POT ")]
                    if len(blocks) != 2 or not all(math.isfinite(float(val)) for val in blocks):
                        raise ValueError((c, blocks))
                    writer.writerow([*c, ",".join(blocks), "0", "0"])
                else:
                    blocks, = [s.split()[1:] for s in proc.stdout.splitlines() if s.startswith("TVN ")]
                    expected_len = 8 + 3 * nc + nc * nc
                    if len(blocks) != expected_len or not all(math.isfinite(float(val)) for val in blocks):
                        raise ValueError((c, len(blocks), expected_len, blocks))
                    writer.writerow([*c, ",".join(blocks[:-2]), blocks[-2], blocks[-1]])
        env["SAFTVRMIE_FIXTURE"] = str(live)
        subprocess.run(["cargo", "test", "--test", "saft_vrmie", "exact_source_tvn"], cwd=ROOT, env=env, check=True)
        if args.write_fixture:
            FIXTURE.write_bytes(live.read_bytes())
            metadata = {"cases": len(all_cases), "fixture_sha256": sha(FIXTURE),
                        "probe_sha256": sha(PROBE), "build": json.loads((ref / "metadata.json").read_text())}
            FIXTURE.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"{len(all_cases)} live Mie cases passed")


if __name__ == "__main__":
    main()
