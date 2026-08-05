#!/usr/bin/env python3
"""Attempt Magma V2.29 saturation of the published rank-29 subgroup through 4095.

The program is submitted to the official free Magma calculator. It is
self-contained and outputs only compact evidence so the 20,000-character cap is
irrelevant. A warning, timeout, Magma error, or missing sentinel is a failure.
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


def program() -> str:
    curve = json.loads((ROOT / "curve.json").read_text())
    pts = json.loads((ROOT / "points.json").read_text())["points"]
    a = curve["a_invariants"]
    points = ",\n  ".join(f"E![Q!({p['x']}),Q!({p['y']}),Q!1]" for p in pts)
    return f'''SetColumns(0);
SetSeed(20260806);
Q := Rationals();
E := EllipticCurve([Q!({a[0]}),Q!({a[1]}),Q!({a[2]}),Q!({a[3]}),Q!({a[4]})]);
Pts := [
  {points}
];
assert #Pts eq 29;
assert &and[P in E : P in Pts];
ME, phi, psi := MinimalModel(E);
assert IsMinimalModel(ME);
MPts := [phi(P) : P in Pts];
assert IsLinearlyIndependent(MPts);
oldH := HeightPairingMatrix(MPts : Precision := 100);
oldR := Determinant(oldH);
S := Saturation(MPts, 4095 : TorsionFree := true, Check := false);
assert #S eq 29;
assert IsLinearlyIndependent(S);
newH := HeightPairingMatrix(S : Precision := 100);
newR := Determinant(newH);
ratio := oldR/newR;
nearest := Round(Sqrt(ratio));
print "R30_MAGMA_VERSION", GetVersion();
print "R30_SATURATION_BOUND", 4095;
print "R30_INPUT_GENERATORS", #MPts;
print "R30_SATURATED_GENERATORS", #S;
print "R30_OLD_REGULATOR", oldR;
print "R30_NEW_REGULATOR", newR;
print "R30_REGULATOR_RATIO", ratio;
print "R30_INDEX_FROM_REGULATOR", nearest;
assert Abs(ratio - nearest^2) lt 10^(-50);
assert nearest eq 1;
print "R30_SATURATION_4095_PASS";
'''


def main() -> int:
    code = program()
    (RESULTS / "magma_saturation_4095.m").write_text(code, encoding="utf-8")
    data = urllib.parse.urlencode({"input": code}).encode()
    request = urllib.request.Request(
        SERVER,
        data=data,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/xml,text/xml,text/html",
            "Referer": REFERER,
            "User-Agent": "rank30-reproducibility-pipeline/saturation-1.0",
        },
        method="POST",
    )
    started = time.time()
    raw = b""
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read()
            http_status = response.status
        (RESULTS / "magma_saturation_4095.xml").write_bytes(raw)
        root = ET.fromstring(raw)
        warning = " ".join("".join(w.itertext()).strip() for w in root.findall(".//warning")).strip()
        lines = ["".join(line.itertext()) for line in root.findall(".//results/line")]
        output = "\n".join(lines)
        (RESULTS / "magma_saturation_4095.txt").write_text(output + "\n", encoding="utf-8")
        error_markers = ("Runtime error", "User error", "Internal error", "Assertion failed", "Identifier '", "Syntax error")
        has_error = any(marker in output for marker in error_markers)
        timed_out = "time limit" in warning.lower() or "timed out" in warning.lower()
        passed = "R30_SATURATION_4095_PASS" in output and not warning and not has_error
        status = "pass" if passed else ("timeout" if timed_out else "fail")
        result = {
            "service": SERVER,
            "status": status,
            "http_status": http_status,
            "elapsed_seconds": time.time() - started,
            "warning": warning,
            "has_error_marker": has_error,
            "sentinel": "R30_SATURATION_4095_PASS",
            "output": output,
            "program_file": "magma_saturation_4095.m",
            "xml_file": "magma_saturation_4095.xml",
        }
    except Exception:
        if raw:
            (RESULTS / "magma_saturation_4095.xml").write_bytes(raw)
        result = {
            "service": SERVER,
            "status": "exception",
            "elapsed_seconds": time.time() - started,
            "traceback": traceback.format_exc(),
        }
    (RESULTS / "magma_saturation_4095.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
