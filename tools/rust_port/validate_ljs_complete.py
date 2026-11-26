"""Validate complete LJS/WCA/UV/UF and cavity Rust parity against live Fortran."""
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
FIXTURE = ROOT / "tests/fixtures/ljs_complete.csv"
PROBE = Path(__file__).with_name("ljs_probe.f90")


def cases():
    # (model, variant, kind, phase_flag, T, v_or_p, n1, sigma, eps)
    def c(model, variant, kind, t, v_or_p, n1=1.2, sigma=0.0, eps=0.0, phase=0):
        return [model, variant, kind, phase, t, v_or_p, n1, sigma, eps]

    # 1. Standalone cavity integral TVN (below eta_c=0.52 and above eta_c=0.52)
    yield c("CAVITY", "default", 6, 120.0, 9.0e-5, 1.2, 3.42e-10, 124.0)
    yield c("CAVITY", "default", 6, 300.0, 1.5e-4, 1.0, 3.42e-10, 124.0)
    yield c("CAVITY", "default", 6, 80.0, 4.3e-5, 1.2, 3.42e-10, 124.0)  # eta > 0.52 extrapolation
    yield c("CAVITY", "default", 6, 45.0, 5.0e-5, 1.2, 2.801e-10, 33.921)

    # 2. Binary calcZetaX_vdW_no_segments
    yield c("ZETAX2", "default", 6, 150.0, 1.2e-4, 1.3, 3.42e-10, 124.0)

    # 3. Standalone UF kernels (LJS-UF and LJ-UF)
    yield c("UF-KERNEL", "ljs", 7, 2.0, 0.25, 1.0, 3.42e-10, 124.0)
    yield c("UF-KERNEL", "lj", 7, 1.5, 0.65, 1.0, 3.42e-10, 124.0)
    yield c("UF-KERNEL", "ljs", 8, 45.0, 6.0e-5, 1.2, 2.801e-10, 33.921)
    yield c("UF-KERNEL", "lj", 8, 120.0, 8.5e-5, 1.3, 3.42e-10, 124.0)

    # 4. LJS-BH: full 10-block TVN, options, custom params, DFT hooks, TP states, model hooks
    yield c("LJS-BH", "default", 0, 100.0, 1.0e-4, 1.2)
    yield c("LJS-BH", "default", 0, 220.0, 8.0e-4, 1.5)
    yield c("LJS-BH", "lafitte", 0, 110.0, 9.0e-5, 1.2)
    yield c("LJS-BH", "no-chi", 0, 110.0, 9.0e-5, 1.2)
    yield c("LJS-BH", "bh-a1-only", 0, 120.0, 1.1e-4, 1.2)
    yield c("LJS-BH", "bh-a2-only", 0, 120.0, 1.1e-4, 1.2)
    yield c("LJS-BH", "bh-a3-only", 0, 120.0, 1.1e-4, 1.2)
    yield c("LJS-BH", "bh-hs-only", 0, 120.0, 1.1e-4, 1.2)
    yield c("LJS-BH", "kolafa", 0, 120.0, 1.1e-4, 1.2)
    yield c("LJS-BH", "default", 0, 150.0, 1.2e-4, 1.1, 3.60e-10, 135.0)
    yield c("LJS-BH", "default", 1, 120.0, 1.0e-4, 1.2)  # calc_ljs_dispersion
    yield c("LJS-BH", "default", 2, 120.0, 1.0e-4, 1.2)  # calc_ljs_hard_sphere
    yield c("LJS-BH", "default", 4, 85.0, 2.0e6, 1.0, phase=1)   # liquid TP
    yield c("LJS-BH", "default", 4, 130.0, 1.0e6, 1.0, phase=2)  # vapor/supercritical TP
    yield c("LJS-BH", "lafitte", 4, 90.0, 3.0e6, 1.0, phase=1)   # Lafitte liquid TP
    yield c("LJS-BH", "default", 5, 120.0, 1.15, 1.0)            # hooks + uncapped potential
    yield c("LJS-BH", "default", 5, 120.0, 0.5, 1.0)             # hooks + capped potential

    # 5. LJS-WCA: full 10-block TVN, all options, eta > 0.52, Kolafa, custom params, DFT hooks, TP
    yield c("LJS-WCA", "default", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-WCA", "default", 0, 220.0, 7.0e-4, 1.5)
    yield c("LJS-WCA", "default", 0, 80.0, 4.3e-5, 1.2)          # eta > 0.52
    yield c("LJS-WCA", "no-cavity", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-WCA", "no-hs", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-WCA", "kolafa", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-WCA", "wca-a1-only", 0, 120.0, 1.0e-4, 1.2)
    yield c("LJS-WCA", "wca-a2-only", 0, 120.0, 1.0e-4, 1.2)
    yield c("LJS-WCA", "wca-a3-only", 0, 120.0, 1.0e-4, 1.2)
    yield c("LJS-WCA", "wca-a4-only", 0, 120.0, 1.0e-4, 1.2)
    yield c("LJS-WCA", "default", 0, 95.0, 8.5e-5, 1.3, 3.55e-10, 118.0)
    yield c("LJS-WCA", "default", 1, 110.0, 9.5e-5, 1.2)         # calc_ljs_dispersion
    yield c("LJS-WCA", "default", 2, 110.0, 9.5e-5, 1.2)         # calc_ljs_hard_sphere
    yield c("LJS-WCA", "default", 3, 110.0, 9.5e-5, 1.2)         # calc_wca_soft_repulsion
    yield c("LJS-WCA", "default", 4, 85.0, 2.0e6, 1.0, phase=1)  # liquid TP
    yield c("LJS-WCA", "default", 4, 130.0, 1.0e6, 1.0, phase=2) # vapor TP
    yield c("LJS-WCA", "default", 4, 250.0, 1.5e7, 1.0, phase=2) # supercritical TP
    yield c("LJS-WCA", "default", 5, 120.0, 1.35, 1.0)           # hooks + spline tail potential

    # 6. LJS-UV: full 10-block TVN, temp-u, no-virial, no-cavity, Kolafa, custom, DFT hooks, TP
    yield c("LJS-UV", "default", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UV", "default", 0, 240.0, 6.0e-4, 1.4)
    yield c("LJS-UV", "default", 0, 80.0, 4.3e-5, 1.2)           # eta > 0.52
    yield c("LJS-UV", "temp-u", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UV", "no-virial", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UV", "no-cavity", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UV", "kolafa", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UV", "default", 0, 140.0, 1.1e-4, 1.1, 3.65e-10, 130.0)
    yield c("LJS-UV", "default", 1, 110.0, 9.5e-5, 1.2)          # calc_ljs_dispersion
    yield c("LJS-UV", "default", 2, 110.0, 9.5e-5, 1.2)          # calc_ljs_hard_sphere
    yield c("LJS-UV", "default", 3, 110.0, 9.5e-5, 1.2)          # calc_wca_soft_repulsion
    yield c("LJS-UV", "default", 4, 85.0, 2.0e6, 1.0, phase=1)   # liquid TP
    yield c("LJS-UV", "default", 4, 130.0, 1.0e6, 1.0, phase=2)  # vapor TP
    yield c("LJS-UV", "temp-u", 4, 130.0, 1.0e6, 1.0, phase=2)   # temp-u TP
    yield c("LJS-UV", "default", 5, 120.0, 1.80, 1.0)            # hooks + beyond cutoff

    # 7. LJS-UF and LJ-UF: full 10-block TVN, no-cavity, Kolafa, DFT hooks, TP
    yield c("LJS-UF", "default", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UF", "no-cavity", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UF", "kolafa", 0, 100.0, 9.0e-5, 1.2)
    yield c("LJS-UF", "default", 1, 110.0, 9.5e-5, 1.2)
    yield c("LJS-UF", "default", 2, 110.0, 9.5e-5, 1.2)
    yield c("LJS-UF", "default", 3, 110.0, 9.5e-5, 1.2)
    yield c("LJS-UF", "default", 4, 85.0, 2.0e6, 1.0, phase=1)
    yield c("LJS-UF", "default", 4, 130.0, 1.0e6, 1.0, phase=2)
    yield c("LJ-UF", "default", 0, 120.0, 9.0e-5, 1.2)
    yield c("LJ-UF", "no-cavity", 0, 120.0, 9.0e-5, 1.2)
    yield c("LJ-UF", "default", 1, 120.0, 9.0e-5, 1.2)
    yield c("LJ-UF", "default", 4, 105.0, 3.0e6, 1.0, phase=1)
    yield c("LJ-UF", "default", 4, 160.0, 1.5e6, 1.0, phase=2)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def verify_upstream_defects(exe: Path, env: dict):
    """Independently verify confirmed upstream Fortran defects live."""
    # 1. ljx_ux_eos_constructor reversed INDEX bug
    stdin = "DEFECTS\nuf-ctor\n8 0\n100.0 0.0 1.0 0.0 0.0\n"
    proc = subprocess.run([exe], input=stdin, text=True, env=env, capture_output=True, check=True)
    line, = [s for s in proc.stdout.splitlines() if s.startswith("DEFECT_UF_CTOR ")]
    flags = [int(x) for x in line.split()[1:]]
    assert flags == [1, 0, 1, 0], f"Unexpected ljx_ux_eos_constructor flags: {flags}"

    # 2. ljs_bh_set_pure_params stale sigma_eff bug
    stdin = "DEFECTS\nbh-set-pure\n8 0\n100.0 0.0 1.0 3.7e-10 140.0\n"
    proc = subprocess.run([exe], input=stdin, text=True, env=env, capture_output=True, check=True)
    line, = [s for s in proc.stdout.splitlines() if s.startswith("DEFECT_BH_SETPARAM ")]
    vals = [float(x) for x in line.split()[1:]]
    assert abs(vals[0] - 3.7e-10) < 1e-22, vals
    assert abs(vals[1] - 3.42e-10) < 1e-22, f"Expected stale sigma_eff=3.42e-10, got {vals[1]}"

    # 3. calc_cavity_integral_LJ_Fres_test zero F_e/F_ee/F_Te bug
    stdin = "DEFECTS\ncavity-test\n8 0\n120.0 0.35 1.0 3.42e-10 124.0\n"
    proc = subprocess.run([exe], input=stdin, text=True, env=env, capture_output=True, check=True)
    line, = [s for s in proc.stdout.splitlines() if s.startswith("DEFECT_CAV_TEST ")]
    vals = [float(x) for x in line.split()[1:]]
    assert vals[0] != 0.0 and vals[2] == 0.0 and vals[4] == 0.0 and vals[5] == 0.0, vals

    # 4. init_ljs('UF') aborts in assocIndices_bookkeeping because eosLJS_UF is missing
    stdin = "DEFECTS\nljs-uf-init\n8 0\n100.0 0.0 1.0 0.0 0.0\n"
    proc = subprocess.run([exe], input=stdin, text=True, env=env, capture_output=True)
    assert proc.returncode != 0 and "At least one CPA-component must self-associate" in proc.stdout, (
        proc.returncode,
        proc.stdout,
        proc.stderr,
    )
    print("Confirmed all 4 live upstream Fortran defects (uf-ctor, bh-set-pure, cavity-test, ljs-uf-init)")


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
        exe = Path(tmp) / ("ljs_probe.exe" if os.name == "nt" else "ljs_probe")
        subprocess.run(
            [compiler, "-fdefault-real-8", "-fdefault-double-8", "-I", ref, PROBE, lib, *flags, "-o", exe],
            check=True,
            env=env,
        )
        verify_upstream_defects(exe, env)
        live = Path(tmp) / "live.csv"
        with live.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter=";", lineterminator="\n")
            writer.writerow("model variant kind phase T v_or_p n1 sigma eps tag values".split())
            for row in all_cases:
                model, variant, kind, phase, t, v_or_p, n1, sigma, eps = row
                stdin = f"{model}\n{variant}\n{kind} {phase}\n{t!r} {v_or_p!r} {n1!r} {sigma!r} {eps!r}\n"
                proc = subprocess.run([exe], input=stdin, text=True, env=env, capture_output=True, timeout=30)
                if proc.returncode:
                    raise RuntimeError((row, proc.stdout, proc.stderr))
                matched = [
                    s.split()
                    for s in proc.stdout.splitlines()
                    if s.startswith(("TVN ", "TP ", "HOOKS ", "ZETAX2 ", "UF6 "))
                ]
                if len(matched) != 1:
                    raise RuntimeError((row, proc.stdout, proc.stderr))
                tag = matched[0][0]
                values = matched[0][1:]
                expected = {"TVN": 10, "TP": 20, "HOOKS": 10, "ZETAX2": 25, "UF6": 6}[tag]
                if len(values) != expected or not all(math.isfinite(float(v)) for v in values):
                    raise ValueError((row, tag, len(values), expected, values))
                writer.writerow([*row, tag, ",".join(values)])
        env["LJS_COMPLETE_FIXTURE"] = str(live)
        subprocess.run(
            ["cargo", "test", "--test", "ljs_complete", "exact_source_ljs_complete"],
            cwd=ROOT,
            env=env,
            check=True,
        )
        if args.write_fixture:
            FIXTURE.write_bytes(live.read_bytes())
            metadata = {
                "cases": len(all_cases),
                "fixture_sha256": sha(FIXTURE),
                "probe_sha256": sha(PROBE),
                "build": json.loads((ref / "metadata.json").read_text(encoding="utf-8")),
            }
            FIXTURE.with_suffix(".metadata.json").write_text(
                json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
            )
    print(f"{len(all_cases)} live LJS/WCA/UV/UF/cavity cases passed")


if __name__ == "__main__":
    main()
