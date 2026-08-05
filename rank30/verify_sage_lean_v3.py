#!/usr/bin/env sage -python
"""Lean, fail-closed SageMath verification of the rank-29 baseline.

Unlike the earlier attempt, this script never asks Sage to rediscover a global
minimal model before doing the independent arithmetic. Magma has separately
certified that the published integral model is already minimal. The modes are
split so that a slow conductor or saturation computation cannot hide the exact
point, finite-reduction, or height evidence.

Inputs are only ``curve.json`` and ``points.json``. Search databases are never
read. Every output is JSON plus a human-readable log under ``results/``.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import subprocess
import time
import traceback
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, QQ, RealField, proof

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
START = time.monotonic()

LOCAL_PRIMES = [19, 23, 29, 37, 47, 53, 59, 73, 79, 83, 97, 101, 103, 107, 109, 127, 131, 151, 157, 163, 173, 179]
TORSION_PRIMES = [19, 23, 67]


def sage_version() -> str:
    try:
        return subprocess.check_output(["sage", "--version"], text=True).strip()
    except Exception:
        return "unknown"


def event(label: str, **data) -> None:
    payload = {"event": label, "elapsed_seconds": round(time.monotonic() - START, 3), **data}
    print("R30_SAGE_LEAN_EVENT " + json.dumps(payload, sort_keys=True), flush=True)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs():
    curve_doc = json.loads((ROOT / "curve.json").read_text(encoding="utf-8"))
    points_doc = json.loads((ROOT / "points.json").read_text(encoding="utf-8"))
    ainvs = [QQ(a) for a in curve_doc["a_invariants"]]
    E = EllipticCurve(QQ, ainvs)
    expected_delta = QQ(curve_doc["discriminant"]["value"])
    assert E.discriminant() == expected_delta
    assert E.discriminant() != 0

    a1, a2, a3, a4, a6 = ainvs
    points = []
    seen = set()
    for rec in points_doc["points"]:
        x = QQ(rec["x"])
        y = QQ(rec["y"])
        assert y * y + a1 * x * y + a3 * y == x**3 + a2 * x**2 + a4 * x + a6
        P = E([x, y, 1])
        assert P in E
        key = (x, y)
        assert key not in seen
        seen.add(key)
        points.append(P)
    assert len(points) == points_doc["count"] == 29
    event("inputs_verified", points=29, discriminant_digits=len(str(abs(int(expected_delta)))))
    return curve_doc, points_doc, E, points


def base_result(curve_doc, E):
    return {
        "system": sage_version(),
        "proof_flags": "proof.all(True)",
        "input_sha256": {
            "curve.json": sha256(ROOT / "curve.json"),
            "points.json": sha256(ROOT / "points.json"),
            "verify_sage_lean_v3.py": sha256(Path(__file__)),
        },
        "curve": {
            "a_invariants": [str(a) for a in E.a_invariants()],
            "discriminant": str(E.discriminant()),
            "nonsingular": bool(E.discriminant() != 0),
            "minimality_note": "This lean Sage run does not invoke global_minimal_model; the same published model is independently certified minimal by Magma V2.29-9.",
        },
    }


def finite_quotient_rows(E, points):
    rows = []
    local_records = []
    for p in LOCAL_PRIMES:
        F = GF(p)
        Ep = E.change_ring(F)
        all_points = list(Ep.points())
        doubles = {2 * R for R in all_points}
        qord = len(all_points) // len(doubles)
        assert qord in (1, 2, 4)
        reductions = [Ep([F(P[0]), F(P[1]), F(1)]) for P in points]
        before = len(rows)
        if qord == 2:
            rows.append([0 if P in doubles else 1 for P in reductions])
        elif qord == 4:
            R1 = next(R for R in all_points if R not in doubles)
            coset1 = {R1 + D for D in doubles}
            R2 = next(R for R in all_points if R not in doubles and R not in coset1)
            coset2 = {R2 + D for D in doubles}
            row1, row2 = [], []
            for P in reductions:
                if P in doubles:
                    bits = (0, 0)
                elif P in coset1:
                    bits = (1, 0)
                elif P in coset2:
                    bits = (0, 1)
                else:
                    bits = (1, 1)
                row1.append(bits[0])
                row2.append(bits[1])
            rows.extend([row1, row2])
        local_records.append({"p": p, "group_order": len(all_points), "quotient_order": qord, "binary_rows_added": len(rows) - before})
    M = Matrix(GF(2), rows)
    return M, local_records


def mode_quick() -> int:
    curve_doc, points_doc, E, points = load_inputs()
    reduction_orders = [int(E.change_ring(GF(p)).cardinality()) for p in TORSION_PRIMES]
    torsion_gcd = math.gcd(*reduction_orders)
    assert reduction_orders == [28, 32, 83]
    assert torsion_gcd == 1
    event("torsion_certified_by_good_reduction", orders=reduction_orders)

    M, local_records = finite_quotient_rows(E, points)
    local_rank = int(M.rank())
    assert M.ncols() == 29
    assert local_rank == 29
    event("local_mod2_independence_certified", rows=M.nrows(), rank=local_rank)

    result = base_result(curve_doc, E)
    result.update({
        "status": "pass",
        "mode": "quick",
        "points_verified": len(points),
        "torsion_certificate": {"good_primes": TORSION_PRIMES, "reduction_orders": reduction_orders, "gcd": torsion_gcd, "conclusion": "E(Q)_tors is trivial"},
        "independence_certificate": {"good_primes": LOCAL_PRIMES, "local_records": local_records, "binary_matrix_rows": M.nrows(), "binary_matrix_columns": M.ncols(), "rank_over_F2": local_rank, "conclusion": "the 29 listed points are Z-linearly independent"},
        "subgroup_rank": 29,
        "unconditional_conclusion": "rank E(Q) >= 29",
    })
    (RESULTS / "sage_lean_v3_quick.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0


def matrix_as_strings(H, digits: int):
    RF = RealField(max(H.base_ring().precision(), 4 * digits + 32))
    return [[str(RF(H[i, j]).n(digits=digits)) for j in range(H.ncols())] for i in range(H.nrows())]


def mode_height(precision: int) -> int:
    curve_doc, points_doc, E, points = load_inputs()
    event("height_start", precision_bits=precision)
    t0 = time.monotonic()
    H = E.height_pairing_matrix(points, precision=precision)
    elapsed = time.monotonic() - t0
    event("height_matrix_computed", precision_bits=precision, seconds=round(elapsed, 3))
    assert H.nrows() == H.ncols() == 29
    assert H == H.transpose()
    det = H.det()
    assert det != 0
    RF = RealField(precision)
    eigenvalues = [RF(v) for v in H.change_ring(RF).eigenvalues()]
    min_eigenvalue = min(eigenvalues)
    max_eigenvalue = max(eigenvalues)
    assert min_eigenvalue > 0
    event("height_positive", determinant=str(det), minimum_eigenvalue=str(min_eigenvalue))

    digits = max(35, int(precision * 0.25))
    result = base_result(curve_doc, E)
    result.update({
        "status": "pass",
        "mode": "height",
        "precision_bits": precision,
        "normalised": True,
        "matrix_dimension": 29,
        "determinant": str(det),
        "minimum_eigenvalue": str(min_eigenvalue),
        "maximum_eigenvalue": str(max_eigenvalue),
        "diagonal": [str(H[i, i]) for i in range(29)],
        "matrix_decimal_digits": digits,
        "matrix": matrix_as_strings(H, digits),
        "rigour_note": "Sage's canonical-height matrix is numerical. Z-independence is not deduced from it here; it is proved exactly by the finite-reduction certificate.",
    })
    json_path = RESULTS / f"sage_lean_v3_height_{precision}.json"
    txt_path = RESULTS / f"sage_lean_v3_height_{precision}.txt"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with txt_path.open("w", encoding="utf-8") as fh:
        fh.write(f"# {sage_version()}\n# precision_bits={precision}\n# determinant={det}\n# minimum_eigenvalue={min_eigenvalue}\n")
        for row in result["matrix"]:
            fh.write(" ".join(row) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "matrix"}, indent=2), flush=True)
    return 0


def mode_saturation(min_prime: int, max_prime: int) -> int:
    curve_doc, points_doc, E, points = load_inputs()
    assert 2 <= min_prime <= max_prime
    event("saturation_start", min_prime=min_prime, max_prime=max_prime)
    log_path = RESULTS / f"sage_lean_v3_saturation_{min_prime}_{max_prime}.log"
    t0 = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print(sage_version(), flush=True)
        print(f"input_generators={len(points)}", flush=True)
        print(f"min_prime={min_prime} max_prime={max_prime}", flush=True)
        saturated_points, index, regulator = E.saturation(points, verbose=True, min_prime=min_prime, max_prime=max_prime)
        print(f"returned_generators={len(saturated_points)}", flush=True)
        print(f"index={index}", flush=True)
        print(f"regulator={regulator}", flush=True)
    elapsed = time.monotonic() - t0
    assert int(index) == 1
    event("saturation_pass", min_prime=min_prime, max_prime=max_prime, index=int(index), seconds=round(elapsed, 3))

    result = base_result(curve_doc, E)
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
        "algorithm_note": "Sage E.saturation uses eclib; index 1 certifies p-saturation for all relevant primes in the requested interval, subject to eclib's computed upper bound.",
    })
    (RESULTS / f"sage_lean_v3_saturation_{min_prime}_{max_prime}.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("quick", "height", "saturation"), required=True)
    ap.add_argument("--precision", type=int, default=128)
    ap.add_argument("--min-prime", type=int, default=2)
    ap.add_argument("--max-prime", type=int, default=4095)
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    proof.all(True)
    try:
        if args.mode == "quick":
            code = mode_quick()
        elif args.mode == "height":
            code = mode_height(args.precision)
        else:
            code = mode_saturation(args.min_prime, args.max_prime)
        raise SystemExit(code)
    except Exception:
        failure = {"status": "fail", "system": sage_version(), "arguments": vars(args), "elapsed_seconds": time.monotonic() - START, "traceback": traceback.format_exc()}
        stem = f"sage_lean_v3_{args.mode.replace('-', '_')}_failure.json"
        (RESULTS / stem).write_text(json.dumps(failure, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(failure, indent=2), flush=True)
        raise
