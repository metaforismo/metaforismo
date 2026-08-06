#!/usr/bin/env sage -python
"""Compute exact local data, conductor, root number, and minimality in Sage.

The complete published factorisation of the discriminant is supplied explicitly,
so this job does not spend hours rediscovering it. Sage still proves each listed
factor prime and runs Tate's algorithm independently at every bad prime.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import time
import traceback
from pathlib import Path

from sage.all import EllipticCurve, QQ, ZZ, proof

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
START = time.monotonic()

FACTORS = [
    (2, 19), (3, 7), (5, 7), (7, 4), (11, 5), (13, 3), (17, 4),
    (31, 3), (41, 2), (43, 2), (61, 2), (233, 1), (241, 2), (4139, 1),
    (678146849364709860535420504397393, 1),
    (159788990966780131363155786084695062643236502969, 1),
    (4402149008473369392540402625019227412319473055901, 1),
]


def version():
    try:
        return subprocess.check_output(["sage", "--version"], text=True).strip()
    except Exception:
        return "unknown"


def event(label, **data):
    print("R30_SAGE_LOCAL_EVENT " + json.dumps({"event": label, "elapsed_seconds": round(time.monotonic() - START, 3), **data}, sort_keys=True), flush=True)


def main() -> int:
    curve_doc = json.loads((ROOT / "curve.json").read_text())
    E = EllipticCurve(QQ, [QQ(a) for a in curve_doc["a_invariants"]])
    delta = ZZ(E.discriminant())
    expected = ZZ(curve_doc["discriminant"]["value"])
    assert delta == expected != 0

    reconstructed = -ZZ(1)
    for p, e in FACTORS:
        event("primality_start", p=str(p), digits=len(str(p)))
        # proof.all(True) is set below, so Integer.is_prime() is a proved test.
        assert ZZ(p).is_prime()
        reconstructed *= ZZ(p) ** e
        event("primality_pass", p=str(p))
    assert reconstructed == delta
    event("factorisation_verified", bad_prime_count=len(FACTORS))

    local_rows = []
    conductor = ZZ(1)
    finite_root = ZZ(1)
    globally_minimal = True
    for p, model_disc_val in FACTORS:
        event("local_data_start", p=str(p))
        data = E.local_data(ZZ(p))
        minimal_disc_val = int(data.discriminant_valuation())
        conductor_val = int(data.conductor_valuation())
        local_root = int(E.root_number(ZZ(p)))
        conductor *= ZZ(p) ** conductor_val
        finite_root *= local_root
        minimal_here = minimal_disc_val == model_disc_val
        globally_minimal = globally_minimal and minimal_here
        row = {
            "p": str(p),
            "digits": len(str(p)),
            "model_discriminant_valuation": model_disc_val,
            "minimal_discriminant_valuation": minimal_disc_val,
            "model_is_minimal_at_p": minimal_here,
            "conductor_valuation": conductor_val,
            "kodaira_symbol": str(data.kodaira_symbol()),
            "tamagawa_number": int(data.tamagawa_number()),
            "bad_reduction_type": data.bad_reduction_type(),
            "local_root_number": local_root,
        }
        local_rows.append(row)
        event("local_data_pass", **row)

    global_root = -finite_root
    assert global_root in (-1, 1)
    assert globally_minimal
    event("global_invariants", conductor=str(conductor), root_number=int(global_root), globally_minimal=True)

    result = {
        "status": "pass",
        "system": version(),
        "proof_flags": "proof.all(True)",
        "inputs": {
            "curve_json_sha256": hashlib.sha256((ROOT / "curve.json").read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "discriminant": str(delta),
        "factorisation": [{"p": str(p), "exponent": e} for p, e in FACTORS],
        "factorisation_verified": True,
        "all_factors_proven_prime": True,
        "local_data": local_rows,
        "global_minimality_certificate": {
            "integral_model": True,
            "minimal_discriminant_valuation_matches_at_every_bad_prime": True,
            "conclusion": "the published model is globally minimal",
        },
        "conductor": str(conductor),
        "finite_local_root_product": int(finite_root),
        "archimedean_root_number": -1,
        "global_root_number": int(global_root),
    }
    (RESULTS / "sage_local_data_v1.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    proof.all(True)
    try:
        raise SystemExit(main())
    except Exception:
        failure = {"status": "fail", "system": version(), "elapsed_seconds": time.monotonic() - START, "traceback": traceback.format_exc()}
        (RESULTS / "sage_local_data_v1_failure.json").write_text(json.dumps(failure, indent=2) + "\n")
        print(json.dumps(failure, indent=2), flush=True)
        raise
