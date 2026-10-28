"""Build and collect isolated, one-request original MBWR coexistence observations."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from reference import ROOT
from validate_meos_flash import digest

HELPER = ROOT / "tools/rust_port/mbwr_coexistence_reference.f90"
FLAGS = ["-O3", "-cpp", "-fdefault-real-8", "-fdefault-double-8", "-frecursive",
         "-fopenmp", "-fallow-argument-mismatch", "-ffree-line-length-none"]
CACHE_SCHEMA = 1
DEPENDENCIES = (ROOT / "tools/rust_port/reference.py",
                ROOT / "tools/rust_port/validate_meos_flash.py")
PER_REQUEST_TIMEOUT = 8.


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _metadata(path):
    return json.loads((path / "metadata.json").read_text(encoding="utf-8"))


def _inputs(build):
    manifest = build / "objects.rsp"
    objects = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        item = line.strip().strip('"')
        if item:
            p = Path(item)
            objects.append((p if p.is_absolute() else build / p).resolve())
    paths = objects + [p.resolve() for p in build.glob("*.mod")] + [manifest.resolve()]
    return {str(p): _sha(p) for p in paths}


def _source_hashes(identity):
    hashes = {}
    for name, expected in identity["sources"].items():
        path = ROOT / "src" / (name + ".f90")
        actual = _sha(path)
        if actual != expected:
            raise RuntimeError(f"original source changed: {path}")
        hashes[name] = actual
    return hashes


def verify_reference(signature):
    """Fail if any original source, helper, manifest, object, or module changed."""
    if digest(HELPER) != signature["helper_sha256"]:
        raise RuntimeError("coexistence reference helper changed")
    if digest(Path(__file__)) != signature["oracle_module_sha256"]:
        raise RuntimeError("coexistence oracle module changed")
    if {p.name: digest(p) for p in DEPENDENCIES} != signature["dependency_hashes"]:
        raise RuntimeError("coexistence oracle dependencies changed")
    if _sha(Path(signature["compiler"])) != signature["compiler_sha256"]:
        raise RuntimeError("compiler binary changed")
    if _source_hashes(signature["reference"]) != signature["source_hashes"]:
        raise RuntimeError("original source hashes changed")
    now = _inputs(Path(signature["build_dir"]))
    if now != signature["object_hashes"]:
        raise RuntimeError("original object/module inputs changed")
    if _sha(Path(signature["reference_library"])) != signature["reference_library_sha256"]:
        raise RuntimeError("reference library changed")


def build_reference(compiler: Path, build_dir: Path):
    compiler, build_dir = Path(compiler).resolve(), Path(build_dir).resolve()
    refdir = ROOT / ".reference"
    identity = _metadata(refdir)
    if identity != _metadata(build_dir):
        raise RuntimeError("reference and source build metadata differ")
    sources = _source_hashes(identity)
    library = refdir / ("thermopack.dll" if os.name == "nt" else "libthermopack.so")
    library_sha = _sha(library)
    if library_sha != identity["library_sha256"]:
        raise RuntimeError("reference library hash differs from metadata")
    objects = _inputs(build_dir)
    env = dict(os.environ)
    env["PATH"] = str(compiler.parent) + os.pathsep + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    version = subprocess.run([str(compiler), "--version"], capture_output=True, text=True,
                             encoding="utf-8", errors="replace", env=env, check=True).stdout.splitlines()[0]
    exe = refdir / ("mbwr_coexistence_reference.exe" if os.name == "nt" else "mbwr_coexistence_reference")
    command = [str(compiler), *FLAGS, "-I", str(build_dir), str(HELPER), "@" + str(build_dir / "objects.rsp"), "-o", str(exe)]
    subprocess.run(command, cwd=refdir, env=env, check=True)
    signature = dict(reference=identity, build_dir=str(build_dir), reference_library=str(library),
                     reference_library_sha256=library_sha, source_hashes=sources, object_hashes=objects,
                     compiler=str(compiler), compiler_sha256=_sha(compiler), compiler_version=version,
                     flags=FLAGS, environment={"OMP_NUM_THREADS": "1"}, helper_sha256=digest(HELPER),
                     oracle_module_sha256=digest(Path(__file__)),
                     dependency_hashes={p.name: digest(p) for p in DEPENDENCIES},
                     per_request_timeout=PER_REQUEST_TIMEOUT, executable_sha256=_sha(exe))
    verify_reference(signature)
    return exe, env, signature


def encode_request(q):
    vals = [q["mode"], q["index"], q["component"], q["T"], q["P"],
            q.get("liquid_volume", 0), q.get("vapor_volume", 0), *q["parameters"],
            *q["critical"], int(bool(q.get("refine", 0)))]
    line = " ".join(map(str, vals)) + "\n"
    if q["mode"] == 5:
        roots = list(q.get("roots", []))
        if len(roots) > 64:
            raise ValueError("mode 5 accepts at most 64 roots")
        line += str(len(roots)) + "\n"
        if roots:
            line += " ".join(map(str, roots)) + "\n"
    return line


def collect_one(exe, env, q, timeout=PER_REQUEST_TIMEOUT):
    started = time.monotonic()
    try:
        run = subprocess.run([str(exe)], input=encode_request(q), text=True, encoding="utf-8",
                             errors="replace", capture_output=True, env=env, timeout=timeout)
        return dict(returncode=run.returncode, stdout=run.stdout, stderr=run.stderr,
                    elapsed_seconds=time.monotonic() - started)
    except subprocess.TimeoutExpired as e:
        decode = lambda x: x.decode("utf-8", errors="replace") if isinstance(x, bytes) else (x or "")
        return dict(censored=True, timeout=timeout, returncode=None, stdout=decode(e.stdout),
                    stderr=decode(e.stderr), elapsed_seconds=time.monotonic() - started)


def _atomic_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    os.replace(tmp, path)


def collect_cached(exe, env, signature, requests, cache_path, reuse=False):
    cache_path = Path(cache_path)
    encoded = [encode_request(q) for q in requests]
    identity = {"schema": CACHE_SCHEMA, "signature": signature,
                "request_bytes": encoded, "request_sha256": hashlib.sha256("".join(encoded).encode("utf-8")).hexdigest()}
    if reuse:
        saved = json.loads(cache_path.read_text(encoding="utf-8"))
        prior = saved.get("identity", {})
        def without_exe(x):
            x = dict(x)
            x.pop("executable_sha256", None)
            # Artifact location can differ between otherwise identical cloned worktrees;
            # the reference library content hash remains part of the signature.
            x.pop("reference_library", None)
            return x
        if prior.get("schema") != CACHE_SCHEMA or without_exe(prior.get("signature", {})) != without_exe(signature):
            raise RuntimeError("stale coexistence collection signature")
        if prior.get("request_bytes") != encoded:
            raise RuntimeError("cached requests differ")
        results = saved.get("results")
        if not isinstance(results, list) or len(results) != len(requests):
            raise RuntimeError("invalid coexistence collection cache")
        historic = prior["signature"]
    else:
        results = []
        historic = signature
        try:
            for i, q in enumerate(requests, 1):
                results.append(collect_one(exe, env, q))
                if i % 25 == 0 or i == len(requests):
                    print(f"collected {i}/{len(requests)}", flush=True)
                if i % 50 == 0 or i == len(requests):
                    _atomic_json(cache_path, {"identity": identity, "results": results})
        finally:
            verify_reference(signature)
    if reuse:
        verify_reference(signature)
    return historic, results
