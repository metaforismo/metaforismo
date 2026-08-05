#!/usr/bin/env sage -python
"""Fail-closed SageMath 10.9 baseline verifier, split into auditable modes.

The script reads only curve.json, points.json, and verify_exact.py.  Expensive
operations are deliberately separated so a conductor factorisation cannot hide
or block the canonical-height and saturation evidence.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import subprocess
import sys
import time
import traceback
from pathlib import Path

from sage.all import EllipticCurve, QQ, RealField, proof

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)


def version() -> str:
    try:
        return subprocess.check_output(["sage", "--version"], text=True).strip()
    except Exception:
        return "unknown"


def mark(label: str, **data) -> None:
    payload = {"event": label, "elapsed_seconds": round(time.monotonic() - START, 3), **data}
    print("R30_SAGE_EVENT " + json.dumps(payload, sort_keys=True), flush=True)


def load_curve_points():
    curve_doc = json.loads((ROOT / "curve.json").read_text(encoding="utf-8"))
    points_doc = json.loads((ROOT / "points.json").read_text(encoding="utf-8"))
    ainvs = [QQ(a) for a in curve_doc["a_invariants"]]
    E = EllipticCurve(QQ, ainvs)
    expected_delta = QQ(curve_doc["discriminant"]["value"])
    assert E.discriminant() == expected_delta != 0
    points = []
    a1, a2, a3, a4, a6 = ainvs
    for i, rec in enumerate(points_doc["points"], start=1):
        x, y = QQ(rec["x"]), QQ(rec["y"])
        assert y * y + a1 * x * y + a3 * y == x**3 + a2 * x**2 + a4 * x + a6
        P = E([x, y, 1])
        assert P in E
        points.append(P)
    assert len(points) == 29
    mark("inputs_verified", points=len(points), discriminant_digits=len(str(abs(int(expected_delta)))))

    Emin = E.global_minimal_model()
    assert Emin == Emin.global_minimal_model()
    if E == Emin:
        points_min = points
        iso = "identity"
    else:
        phi = E.isomorphism_to(Emin)
        points_min = [phi(P) for P in points]
        iso = str(phi)
    assert all(P in Emin for P in points_min)
    mark("minimal_model_ready", same_model=bool(E == Emin))
    return curve_doc, points_doc, E, Emin, points_min, iso


def base_result(curve_doc, points_doc, E, Emin, iso):
    return {
        "system": version(),
        "input_sha256": {
            "curve.json": hashlib.sha256((ROOT / "curve.json").read_bytes()).hexdigest(),
            "points.json": hashlib.sha256((ROOT / "points.json").read_bytes()).hexdigest(),
            "verify_exact.py": hashlib.sha256((ROOT / "verify_exact.py").read_bytes()).hexdigest(),
            "verify_sage_baseline_v2.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "published_model": {
            "a_invariants": [str(a) for a in E.a_invariants()],
            "discriminant": str(E.discriminant()),
        },
        "minimal_model": {
            "a_invariants": [str(a) for a in Emin.a_invariants()],
            "discriminant": str(Emin.discriminant()),
            "isomorphism_from_published_model": iso,
            "idempotent_check": True,
        },
        "points_verified": len(points_doc["points"]),
    }


def mode_quick() -> int:
    curve_doc, points_doc, E, Emin, points, iso = load_curve_points()
    torsion = Emin.torsion_subgroup()
    assert int(torsion.order()) == 1
    mark("torsion_verified", order=int(torsion.order()))

    exact = subprocess.run(
        [sys.executable, str(ROOT / "verify_exact.py")],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    (RESULTS / "sage_v2_exact_subprocess.log").write_text(
        exact.stdout + exact.stderr, encoding="utf-8"
    )
    assert exact.returncode == 0
    mark("exact_local_certificate_verified")

    result = base_result(curve_doc, points_doc, E, Emin, iso)
    result.update({
        "status": "pass",
        "mode": "quick",
        "torsion_order": int(torsion.order()),
        "torsion_invariants": [int(v) for v in torsion.invariants()],
        "exact_finite_reduction_independence": "pass",
        "unconditional_conclusion": "the 29 listed rational points are Z-linearly independent",
    })
    path = RESULTS / "sage_v2_quick.json"
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0


def decimal_matrix(M, digits: int):
    RF = RealField(max(128, int(digits * 3.5) + 32))
    return [[str(RF(M[i, j]).n(digits=digits)) for j in range(M.ncols())]
            for i in range(M.nrows())]


def mode_height(precision: int) -> int:
    curve_doc, points_doc, E, Emin, points, iso = load_curve_points()
    mark("height_matrix_start", precision_bits=precision)
    t0 = time.monotonic()
    H = Emin.height_pairing_matrix(points, precision=precision)
    mark("height_matrix_computed", precision_bits=precision,
         seconds=round(time.monotonic() - t0, 3))
    assert H.nrows() == H.ncols() == 29
    assert H == H.transpose()
    determinant = H.det()
    assert determinant != 0
    RF = RealField(precision)
    eigenvalues = [RF(v) for v in H.change_ring(RF).eigenvalues()]
    minimum_eigenvalue = min(eigenvalues)
    assert minimum_eigenvalue > 0
    mark("height_matrix_positive", determinant=str(determinant),
         minimum_eigenvalue=str(minimum_eigenvalue))

    digits = max(40, int(precision * 0.25))
    matrix = decimal_matrix(H, digits)
    result = base_result(curve_doc, points_doc, E, Emin, iso)
    result.update({
        "status": "pass",
        "mode": "height",
        "precision_bits": precision,
        "normalisation": "SageMath canonical height pairing (normalised=True; over Q this is the BSD normalisation)",
        "determinant": str(determinant),
        "minimum_eigenvalue": str(minimum_eigenvalue),
        "maximum_eigenvalue": str(max(eigenvalues)),
        "diagonal": [str(H[i, i]) for i in range(29)],
        "matrix_decimal_digits": digits,
        "matrix": matrix,
        "independence_note": "positive-definiteness is numerically checked; exact Z-independence is certified separately by finite reductions",
    })
    json_path = RESULTS / f"sage_v2_height_{precision}.json"
    txt_path = RESULTS / f"sage_v2_height_{precision}.txt"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with txt_path.open("w", encoding="utf-8") as fh:
        fh.write(f"# system={version()}\n")
        fh.write(f"# precision_bits={precision}\n")
        fh.write(f"# determinant={determinant}\n")
        fh.write(f"# minimum_eigenvalue={minimum_eigenvalue}\n")
        for row in matrix:
            fh.write(" ".join(row) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "matrix"}, indent=2), flush=True)
    return 0


def mode_conductor_root() -> int:
    curve_doc, points_doc, E, Emin, points, iso = load_curve_points()
    mark("conductor_start")
    t0 = time.monotonic()
    conductor = int(Emin.conductor())
    mark("conductor_done", seconds=round(time.monotonic() - t0, 3), conductor=str(conductor))
    t1 = time.monotonic()
    root_number = int(Emin.root_number())
    mark("root_number_done", seconds=round(time.monotonic() - t1, 3), root_number=root_number)
    result = base_result(curve_doc, points_doc, E, Emin, iso)
    result.update({
        "status": "pass",
        "mode": "conductor_root",
        "conductor": str(conductor),
        "root_number": root_number,
    })
    (RESULTS / "sage_v2_conductor_root.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2), flush=True)
    return 0


def mode_saturation(min_prime: int, max_prime: int) -> int:
    curve_doc, points_doc, E, Emin, points, iso = load_curve_points()
    assert 2 <= min_prime <= max_prime
    mark("saturation_start", min_prime=min_prime, max_prime=max_prime)
    log_path = RESULTS / f"sage_v2_saturation_{min_prime}_{max_prime}.log"
    t0 = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print(version(), flush=True)
        print(f"min_prime={min_prime} max_prime={max_prime}", flush=True)
        saturated_points, index, regulator = Emin.saturation(
            points,
            min_prime=min_prime,
            max_prime=max_prime,
            verbose=True,
        )
        print(f"returned_generators={len(saturated_points)}", flush=True)
        print(f"index={index}", flush=True)
        print(f"regulator={regulator}", flush=True)
    elapsed = time.monotonic() - t0
    assert int(index) == 1
    mark("saturation_done", min_prime=min_prime, max_prime=max_prime,
         seconds=round(elapsed, 3), index=int(index))
    result = base_result(curve_doc, points_doc, E, Emin, iso)
    result.update({
        "status": "pass",
        "mode": "saturation",
        "min_prime": min_prime,
        "max_prime": max_prime,
        "input_generators": len(points),
        "returned_generators": len(saturated_points),
        "index": int(index),
        "regulator": str(regulator),
        "log_file": log_path.name,
    })
    (RESULTS / f"sage_v2_saturation_{min_prime}_{max_prime}.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2), flush=True)
    return 0


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("quick", "height", "conductor-root", "saturation"), required=True)
    ap.add_argument("--precision", type=int, default=128)
    ap.add_argument("--min-prime", type=int, default=2)
    ap.add_argument("--max-prime", type=int, default=4095)
    return ap.parse_args()


START = time.monotonic()

if __name__ == "__main__":
    args = parse_args()
    proof.all(True)
    try:
        if args.mode == "quick":
            rc = mode_quick()
        elif args.mode == "height":
            rc = mode_height(args.precision)
        elif args.mode == "conductor-root":
            rc = mode_conductor_root()
        else:
            rc = mode_saturation(args.min_prime, args.max_prime)
        raise SystemExit(rc)
    except Exception:
        failure = {
            "status": "fail",
            "system": version(),
            "arguments": vars(args),
            "elapsed_seconds": time.monotonic() - START,
            "traceback": traceback.format_exc(),
        }
        name = args.mode.replace("-", "_")
        (RESULTS / f"sage_v2_{name}_failure.json").write_text(
            json.dumps(failure, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(failure, indent=2), flush=True)
        raise
