#!/usr/bin/env python3
"""Strict independent checks via the official free Magma calculator.

Each Magma program is self-contained. A task passes only when its sentinel is
present and the XML response contains neither a warning nor a Magma error.
The calculator's 60-second cap is recorded as a timeout, never as evidence.
"""
from __future__ import annotations

import json
import time
import traceback
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
SERVER = "http://magma.maths.usyd.edu.au/xml/calculator.xml"
REFERER = "http://magma.maths.usyd.edu.au/calc/"


def load_points_code() -> str:
    points_doc = json.loads((ROOT / "points.json").read_text())
    rows = [f"E![Q!({r['x']}),Q!({r['y']}),Q!1]" for r in points_doc["points"]]
    return "Pts := [\n  " + ",\n  ".join(rows) + "\n];\n"


def prelude() -> str:
    curve = json.loads((ROOT / "curve.json").read_text())
    a = curve["a_invariants"]
    return f'''SetColumns(0);
SetSeed(20260806);
Q := Rationals();
E := EllipticCurve([Q!({a[0]}),Q!({a[1]}),Q!({a[2]}),Q!({a[3]}),Q!({a[4]})]);
EXPECTED_DELTA := Integers()!({curve["discriminant"]["value"]});
assert Integers()!Discriminant(E) eq EXPECTED_DELTA;
assert Discriminant(E) ne 0;
''' + load_points_code() + '''
assert #Pts eq 29;
assert &and[P in E : P in Pts];
'''


def task_model_torsion() -> str:
    return prelude() + r'''
ME, phi, psi := MinimalModel(E);
assert IsMinimalModel(ME);
MPts := [phi(P) : P in Pts];
assert &and[P in ME : P in MPts];
T, tmap := TorsionSubgroup(ME);
assert #T eq 1;
print "R30_MAGMA_VERSION", GetVersion();
print "R30_MINIMAL_AINVS", aInvariants(ME);
print "R30_MINIMAL_DISCRIMINANT", Discriminant(ME);
print "R30_TORSION_ORDER", #T;
print "R30_MODEL_TORSION_PASS";
'''


def task_local_mod2() -> str:
    return prelude() + r'''
primes := [19,23,29,37,47,53,59,73,79,83,97,101,103,107,109,127,131,151,157,163,173,179];
torsion_primes := [19,23,67];
orders := [];
for p in torsion_primes do
    Ep := ChangeRing(E,GF(p));
    Append(~orders,#Ep);
end for;
print "R30_TORSION_REDUCTION_ORDERS", orders;
assert GCD(orders) eq 1;
rows := [];
for p in primes do
    F := GF(p);
    Ep := ChangeRing(E,F);
    allp := Setseq(Points(Ep));
    doubles := {2*R : R in allp};
    qord := #allp div #doubles;
    assert qord in {1,2,4};
    red := [Ep![F!P[1],F!P[2],F!1] : P in Pts];
    if qord eq 2 then
        Append(~rows,Vector(GF(2),[P in doubles select 0 else 1 : P in red]));
    elif qord eq 4 then
        R1 := Rep({R : R in allp | not R in doubles});
        coset1 := {R1 + D : D in doubles};
        R2 := Rep({R : R in allp | not R in doubles and not R in coset1});
        coset2 := {R2 + D : D in doubles};
        row1 := []; row2 := [];
        for P in red do
            if P in doubles then
                Append(~row1,0); Append(~row2,0);
            elif P in coset1 then
                Append(~row1,1); Append(~row2,0);
            elif P in coset2 then
                Append(~row1,0); Append(~row2,1);
            else
                Append(~row1,1); Append(~row2,1);
            end if;
        end for;
        Append(~rows,Vector(GF(2),row1));
        Append(~rows,Vector(GF(2),row2));
    end if;
end for;
M := Matrix(rows);
print "R30_LOCAL_ROWS", Nrows(M);
print "R30_LOCAL_MOD2_RANK", Rank(M);
assert Rank(M) eq 29;
print "R30_EXACT_LOCAL_INDEPENDENCE_PASS";
'''


def task_independence() -> str:
    return prelude() + r'''
ME, phi, psi := MinimalModel(E);
MPts := [phi(P) : P in Pts];
b := IsLinearlyIndependent(MPts);
print "R30_IS_LINEARLY_INDEPENDENT", b;
assert b;
print "R30_DIRECT_INDEPENDENCE_PASS";
'''


def task_height(precision: int) -> str:
    return prelude() + f'''
ME, phi, psi := MinimalModel(E);
MPts := [phi(P) : P in Pts];
H := HeightPairingMatrix(MPts : Precision := {precision});
R := Determinant(H);
print "R30_HEIGHT_PRECISION", {precision};
print "R30_REGULATOR", R;
print "R30_HEIGHT_DIAGONAL", [H[i,i] : i in [1..Nrows(H)]];
print "R30_HEIGHT_MATRIX_BEGIN";
print H;
print "R30_HEIGHT_MATRIX_END";
assert Nrows(H) eq 29 and Ncols(H) eq 29;
assert not IsZero(R);
print "R30_HEIGHT_{precision}_PASS";
'''


def task_conductor() -> str:
    return prelude() + r'''
ME := MinimalModel(E);
print "R30_CONDUCTOR", Conductor(ME);
print "R30_CONDUCTOR_PASS";
'''


def task_root_number() -> str:
    return prelude() + r'''
ME := MinimalModel(E);
print "R30_ROOT_NUMBER", RootNumber(ME);
print "R30_ROOT_NUMBER_PASS";
'''


def submit(name: str, code: str, sentinel: str) -> dict:
    started = time.time()
    data = urllib.parse.urlencode({"input": code}).encode()
    request = urllib.request.Request(
        SERVER,
        data=data,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/xml,text/xml,text/html",
            "Referer": REFERER,
            "User-Agent": "rank30-reproducibility-pipeline/2.0",
        },
        method="POST",
    )
    raw = b""
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read()
            status_code = response.status
        (RESULTS / f"magma_{name}.xml").write_bytes(raw)
        root = ET.fromstring(raw)
        warning = " ".join(
            "".join(w.itertext()).strip() for w in root.findall(".//warning")
        ).strip()
        lines = ["".join(line.itertext()) for line in root.findall(".//results/line")]
        output = "\n".join(lines)
        (RESULTS / f"magma_{name}.txt").write_text(output + "\n", encoding="utf-8")
        error_markers = (
            "Runtime error", "User error", "Internal error", "Assertion failed",
            "Identifier '", "Syntax error", ">>",
        )
        has_error = any(marker in output for marker in error_markers)
        timed_out = "time limit" in warning.lower() or "timed out" in warning.lower()
        passed = sentinel in output and not warning and not has_error
        status = "pass" if passed else ("timeout" if timed_out else "fail")
        return {
            "name": name,
            "status": status,
            "sentinel": sentinel,
            "http_status": status_code,
            "elapsed_seconds": time.time() - started,
            "warning": warning,
            "has_error_marker": has_error,
            "output_file": f"magma_{name}.txt",
            "xml_file": f"magma_{name}.xml",
            "output_tail": output[-6000:],
        }
    except Exception:
        if raw:
            (RESULTS / f"magma_{name}.xml").write_bytes(raw)
        return {
            "name": name,
            "status": "exception",
            "elapsed_seconds": time.time() - started,
            "traceback": traceback.format_exc(),
        }


def main() -> int:
    tasks = [
        ("model_torsion", task_model_torsion(), "R30_MODEL_TORSION_PASS", True),
        ("local_mod2", task_local_mod2(), "R30_EXACT_LOCAL_INDEPENDENCE_PASS", True),
        ("independence", task_independence(), "R30_DIRECT_INDEPENDENCE_PASS", True),
        ("height_80", task_height(80), "R30_HEIGHT_80_PASS", True),
        ("height_160", task_height(160), "R30_HEIGHT_160_PASS", True),
        ("conductor", task_conductor(), "R30_CONDUCTOR_PASS", False),
        ("root_number", task_root_number(), "R30_ROOT_NUMBER_PASS", False),
    ]
    results = []
    required = {}
    for name, code, sentinel, is_required in tasks:
        (RESULTS / f"magma_{name}.m").write_text(code, encoding="utf-8")
        result = submit(name, code, sentinel)
        result["required"] = is_required
        results.append(result)
        if is_required:
            required[name] = result["status"] == "pass"
        print(json.dumps(result, indent=2))
    doc = {
        "service": SERVER,
        "strict_parser": True,
        "calculator_limit_note": "Warning, timeout, or Magma error makes a task fail regardless of sentinels.",
        "tasks": results,
        "required_passes": required,
        "core_pass": all(required.values()),
    }
    (RESULTS / "magma_free_result_v2.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8"
    )
    return 0 if doc["core_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
