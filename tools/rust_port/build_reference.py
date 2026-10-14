"""Build this fork with gfortran, using the upstream source dependency order.

Usage: python tools/rust_port/build_reference.py --compiler /path/to/gfortran
Requires the pinned external/lapack submodule. Artifacts stay in .reference/.
Only LAPACK/BLAS routines actually referenced by ThermoPack are compiled.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def build(compiler):
    compiler = Path(compiler).resolve()
    out = ROOT / ".reference"
    out.mkdir(exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = str(compiler.parent) + os.pathsep + env["PATH"]
    version = subprocess.check_output([str(compiler), "--version"], env=env, text=True)
    if re.search(r"\b15\.\d+", version.splitlines()[0]):
        raise RuntimeError("Upstream warns about gfortran 15; use 14 or another qualified version")
    make = (ROOT / "Makefile.code").read_text()
    block = make.split("OBJE =", 1)[1].split("\n\n", 1)[0]
    names = re.findall(r"\$\(ODIR\)/(\w+)\.o", block)
    flags = ["-O3", "-cpp", "-fdefault-real-8", "-fdefault-double-8",
             "-frecursive", "-fopenmp", "-fallow-argument-mismatch", "-ffree-line-length-none"]
    pic_flags = ["-fPIC"] if os.name != "nt" else []
    flags += pic_flags
    source_hashes = {name: hashlib.sha256((ROOT / "src" / (name + ".f90")).read_bytes()).hexdigest()
                     for name in names}
    lapack_commit = subprocess.check_output(["git", "-C", str(ROOT / "external/lapack"),
                                            "rev-parse", "HEAD"], text=True).strip()
    fingerprint = {"compiler": version.splitlines()[0], "flags": flags,
                   "sources": source_hashes, "lapack_commit": lapack_commit}
    marker = out / "build_fingerprint.json"
    old = None
    if marker.exists():
        old = json.loads(marker.read_text())
    elif (out / "metadata.json").exists():
        metadata = json.loads((out / "metadata.json").read_text())
        old = {key: metadata.get(key) for key in fingerprint}
    if old != fingerprint:
        # One directory, explicit file patterns, no recursive filesystem changes.
        for pattern in ["tp_*.o", "la_*.o", "*.mod", "*.smod"]:
            for artifact in out.glob(pattern):
                if artifact.is_file() and artifact.parent.resolve() == out.resolve():
                    artifact.unlink()
    marker.write_text(json.dumps(fingerprint, indent=2) + "\n")
    objects = []
    for i, name in enumerate(names):
        source = ROOT / "src" / (name + ".f90")
        obj = out / ("tp_" + name + ".o")
        if not obj.exists() or obj.stat().st_mtime < source.stat().st_mtime:
            subprocess.run([str(compiler), *flags, "-J", str(out), "-I", str(out),
                            "-c", str(source), "-o", str(obj)], cwd=out, env=env, check=True)
        objects.append(obj)
        if i % 15 == 0:
            print(f"ThermoPack {i + 1}/{len(names)}: {name}", flush=True)

    lapack = ROOT / "external/lapack"
    candidates = {}
    for directory in [lapack / "SRC", lapack / "BLAS/SRC", lapack / "INSTALL"]:
        for source in sorted(directory.glob("*")):
            if source.suffix.lower() in (".f", ".f90"):
                candidates[source.stem.lower()] = source
    if not candidates:
        raise RuntimeError("Initialize external/lapack first")
    for name in ["la_constants", "la_xisnan"]:
        if name in candidates:
            obj = out / ("la_" + name + ".o")
            if not obj.exists():
                subprocess.run([str(compiler), "-O3", "-cpp", *pic_flags,
                                "-J", str(out), "-I", str(out), "-c", str(candidates[name]),
                                "-o", str(obj)], env=env, cwd=out, check=True)
            objects.append(obj)
    nm = compiler.parent / ("nm.exe" if os.name == "nt" else "nm")
    def symbols(paths):
        defined, undefined = set(), set()
        for offset in range(0, len(paths), 25):
            output = subprocess.check_output([str(nm), "-g", *map(str, paths[offset:offset+25])],
                                             env=env, text=True)
            for line in output.splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[-2] == "U":
                    undefined.add(parts[-1])
                elif len(parts) >= 3 and parts[-2] in ("T", "D", "B", "R"):
                    defined.add(parts[-1])
        return undefined - defined

    compiled = set()
    while True:
        required = sorted({symbol[:-1] for symbol in symbols(objects)
                           if symbol.endswith("_") and symbol[:-1] in candidates} - compiled)
        if not required:
            break
        for name in required:
            obj = out / ("la_" + name + ".o")
            if not obj.exists():
                subprocess.run([str(compiler), "-O3", "-frecursive", "-fallow-argument-mismatch", *pic_flags,
                                "-I", str(out), "-c", str(candidates[name]), "-o", str(obj)],
                               env=env, cwd=out, check=True)
            objects.append(obj)
            compiled.add(name)
        print(f"LAPACK/BLAS: {len(compiled)} routines", flush=True)

    response = out / "objects.rsp"
    response.write_text("\n".join('"' + p.as_posix() + '"' for p in objects))
    library = out / ("thermopack.dll" if os.name == "nt" else "libthermopack.so")
    subprocess.run([str(compiler), "-shared", "-fopenmp", "@" + str(response),
                    "-o", str(library)], cwd=out, env=env, check=True)
    metadata = {"baseline": "d68c794", "compiler": version.splitlines()[0], "flags": flags,
                "lapack_commit": lapack_commit,
                "library_sha256": hashlib.sha256(library.read_bytes()).hexdigest(),
                "sources": source_hashes}
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Built {library}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", required=True)
    build(parser.parse_args().compiler)
