"""Exact-source parity qualification for the Lee-Kesler equation of state."""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from deadline_reference import DeadlineReferenceProcess
from reference import ROOT


def cases():
    grid = [
        # (components, parameter_reference, composition, [(T, P, phase), ...])
        ("C1", "Default", [1.0], [(150.0, 2.0e6, "L"), (250.0, 1.0e6, "V"), (300.0, 5.0e6, "V")]),
        ("C3", "Default", [1.0], [(280.0, 2.0e6, "L"), (350.0, 1.0e6, "V"), (420.0, 6.0e6, "V")]),
        ("CO2", "Default", [1.0], [(260.0, 4.0e6, "L"), (320.0, 2.0e6, "V"), (350.0, 1.0e7, "L")]),
        (
            "C1,C2",
            "Default",
            [0.6, 0.4],
            [(180.0, 3.0e6, "L"), (250.0, 2.0e6, "V"), (350.0, 5.0e6, "V")],
        ),
        (
            "C1,C2",
            "Plocker1982(+)",
            [0.6, 0.4],
            [(200.0, 1.0e6, "L"), (260.0, 2.0e6, "V"), (350.0, 5.0e6, "V")],
        ),
        (
            "C1,C3",
            "Default",
            [0.7, 0.3],
            [(220.0, 4.0e6, "L"), (320.0, 2.0e6, "V"), (380.0, 8.0e6, "V")],
        ),
        (
            "CO2,C1",
            "Default",
            [0.5, 0.5],
            [(240.0, 5.0e6, "L"), (310.0, 4.0e6, "V"), (360.0, 1.0e7, "V")],
        ),
        (
            "C1,C2,C3",
            "Default",
            [0.5, 0.3, 0.2],
            [(200.0, 4.0e6, "L"), (300.0, 2.0e6, "V"), (380.0, 8.0e6, "V")],
        ),
    ]
    for comps, ref, z, states in grid:
        for t, p, phase in states:
            yield comps, ref, z, t, p, phase


def validate(directory):
    good = []
    failed = []
    censored = []
    unqualified = []
    with DeadlineReferenceProcess(directory) as oracle:
        for comps, ref, z, t, p, phase in cases():
            label = f"LEE_KESLER|{comps}|{ref}|{t}|{p}|{','.join(map(repr, z))}|{phase}"
            request = dict(
                mode="LEE_KESLER",
                components=comps,
                parameter_reference=ref,
                feed=z,
                temperature=t,
                pressure=p,
                phase=phase,
            )
            result = oracle.call(request, timeout=10.0)
            record = dict(case=label, request=request, result=result)
            if "reference_timeout" in result:
                censored.append(record)
                continue
            if "values" not in result:
                failed.append(record)
                continue
            if not all(map(math.isfinite, result["values"])):
                unqualified.append(record)
                continue
            good.append(record)

    subprocess.run(
        ["cargo", "build", "--release", "--example", "parity", "-j", "2"],
        cwd=ROOT,
        check=True,
    )
    exe = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    output = subprocess.run(
        [str(exe)],
        input="\n".join(r["case"] for r in good) + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    assert len(output) == len(good)

    rows = []
    differences = []
    worst_rel = 0.0
    rtol = 1e-8
    atol = 1e-10
    for r, out in zip(good, output):
        expected = r["result"]["values"]
        rows.append([r["case"], *[f"{v:.17e}" for v in expected]])
        if out.startswith("ERROR"):
            differences.append(dict(case=r["case"], rust=out))
            continue
        native = list(map(float, out.split(",")))
        assert len(native) == len(expected)
        bad = []
        for idx, (a, b) in enumerate(zip(expected, native)):
            scale = max(abs(a), abs(b), 1.0)
            err = abs(a - b)
            worst_rel = max(worst_rel, err / scale)
            if not math.isfinite(b) or err > atol + rtol * max(abs(a), abs(b)):
                bad.append([idx, a, b])
        if bad:
            differences.append(dict(case=r["case"], quantities=bad))

    path = ROOT / "tests/fixtures/lee_kesler.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "values"])
        writer.writerows(rows)

    meta = json.loads((ROOT / ".reference/metadata.json").read_text(encoding="utf-8"))
    meta.update(
        cases=len(rows),
        generator="tools/rust_port/validate_lee_kesler.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=rtol,
        absolute_tolerance=atol,
        reference_errors=failed,
        censored_calls=censored,
        unqualified_returns=unqualified,
        worst_scaled_error=worst_rel,
    )
    path.with_suffix(".metadata.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    (ROOT / ".reference/lee-kesler-differences.json").write_text(
        json.dumps(differences, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            dict(
                qualified=len(rows),
                failed=len(failed),
                censored=len(censored),
                unqualified=len(unqualified),
                differences=len(differences),
                worst_scaled_error=worst_rel,
            ),
            indent=2,
        )
    )
    if differences:
        print(json.dumps(differences[:8], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    validate(parser.parse_args().compiler_directory)
