#!/usr/bin/env sage -python
"""Directly prove p-saturation of the 29-point subgroup for a prime interval.

This avoids the global height-bound stage in ``E.saturation``.  For every prime
p in the requested interval it calls SageMath's independently implemented
``EllipticCurveSaturator.p_saturation(..., sieve=True)``.  The documented return
value ``False`` is an exact certificate that the given subgroup is p-saturated;
any replacement point, exception, timeout, or missing result fails closed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
import traceback
from pathlib import Path

from sage.all import EllipticCurve, QQ, prime_range, proof
from sage.schemes.elliptic_curves.saturation import EllipticCurveSaturator

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
START = time.monotonic()


def version() -> str:
    try:
        return subprocess.check_output(["sage", "--version"], text=True).strip()
    except Exception:
        return "unknown"


def event(label: str, **data) -> None:
    payload = {"event": label, "elapsed_seconds": round(time.monotonic() - START, 3), **data}
    print("R30_SAGE_PSAT_EVENT " + json.dumps(payload, sort_keys=True, default=str), flush=True)


def load_curve_points():
    curve_doc = json.loads((ROOT / "curve.json").read_text(encoding="utf-8"))
    points_doc = json.loads((ROOT / "points.json").read_text(encoding="utf-8"))
    E = EllipticCurve(QQ, [QQ(a) for a in curve_doc["a_invariants"]])
    assert E.discriminant() == QQ(curve_doc["discriminant"]["value"])
    points = [E([QQ(rec["x"]), QQ(rec["y"]), QQ(1)]) for rec in points_doc["points"]]
    assert len(points) == 29 and all(P in E for P in points)
    return E, points


def main(min_prime: int, max_prime: int) -> int:
    assert 2 <= min_prime <= max_prime
    E, points = load_curve_points()
    primes = [int(p) for p in prime_range(min_prime, max_prime + 1)]
    assert primes
    event("interval_start", min_prime=min_prime, max_prime=max_prime, prime_count=len(primes))
    saturator = EllipticCurveSaturator(E, verbose=False)
    records = []
    for p in primes:
        started = time.monotonic()
        result = saturator.p_saturation(points, p, sieve=True)
        elapsed = time.monotonic() - started
        saturated = result is False
        record = {"p": p, "saturated": saturated, "seconds": round(elapsed, 6)}
        records.append(record)
        event("prime_done", **record)
        if not saturated:
            raise RuntimeError(f"Subgroup is not {p}-saturated: {result!r}")

    output = {
        "status": "pass",
        "system": version(),
        "proof_flags": "proof.all(True)",
        "method": "EllipticCurveSaturator.p_saturation with sieve=True",
        "min_prime": min_prime,
        "max_prime": max_prime,
        "prime_count": len(primes),
        "primes": primes,
        "records": records,
        "all_p_saturated": True,
        "input_point_count": 29,
        "input_sha256": {
            "curve.json": hashlib.sha256((ROOT / "curve.json").read_bytes()).hexdigest(),
            "points.json": hashlib.sha256((ROOT / "points.json").read_bytes()).hexdigest(),
            "script": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
    }
    path = RESULTS / f"sage_direct_p_saturation_{min_prime}_{max_prime}.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-prime", type=int, required=True)
    parser.add_argument("--max-prime", type=int, required=True)
    args = parser.parse_args()
    proof.all(True)
    try:
        raise SystemExit(main(args.min_prime, args.max_prime))
    except Exception:
        failure = {
            "status": "fail",
            "system": version(),
            "arguments": vars(args),
            "elapsed_seconds": time.monotonic() - START,
            "traceback": traceback.format_exc(),
        }
        (RESULTS / f"sage_direct_p_saturation_{args.min_prime}_{args.max_prime}_failure.json").write_text(
            json.dumps(failure, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(failure, indent=2), flush=True)
        raise
