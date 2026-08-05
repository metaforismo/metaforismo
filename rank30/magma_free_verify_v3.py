#!/usr/bin/env python3
"""Fail-closed independent Magma V2.29 verification with bounded output.

The official free Magma calculator caps each program at 60 seconds and each
response at 20,000 characters.  The height matrix is therefore recomputed in
small row blocks.  Every task is self-contained and passes only when its
sentinel is present and the XML response has no warning or Magma error.
"""
from __future__ import annotations

import decimal
import json
import re
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
    points_doc = json.loads((ROOT / "points.json").read_text(encoding="utf-8"))
    rows = [f"E![Q!({r['x']}),Q!({r['y']}),Q!1]" for r in points_doc["points"]]
    return "Pts := [\n  " + ",\n  ".join(rows) + "\n];\n"


def prelude() -> str:
    curve = json.loads((ROOT / "curve.json").read_text(encoding="utf-8"))
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
print "R30_TORSION_REDUCTION_ORDERS", orders;
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


def task_height_summary(precision: int) -> str:
    return prelude() + f'''
ME, phi, psi := MinimalModel(E);
MPts := [phi(P) : P in Pts];
H := HeightPairingMatrix(MPts : Precision := {precision});
R := Determinant(H);
assert Nrows(H) eq 29 and Ncols(H) eq 29;
assert not IsZero(R);
print "R30_HEIGHT_PRECISION", {precision};
print "R30_REGULATOR", R;
print "R30_HEIGHT_DIAGONAL", [H[i,i] : i in [1..29]];
print "R30_HEIGHT_SUMMARY_{precision}_PASS";
'''


def task_height_rows(precision: int, first: int, last: int) -> str:
    return prelude() + f'''
ME, phi, psi := MinimalModel(E);
MPts := [phi(P) : P in Pts];
H := HeightPairingMatrix(MPts : Precision := {precision});
assert Nrows(H) eq 29 and Ncols(H) eq 29;
for i in [{first}..{last}] do
    print "R30_HEIGHT_ROW", i, [H[i,j] : j in [1..29]];
end for;
print "R30_HEIGHT_ROWS_{precision}_{first}_{last}_PASS";
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
            "User-Agent": "rank30-reproducibility-pipeline/3.0",
        },
        method="POST",
    )
    raw = b""
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read()
            status_code = response.status
        (RESULTS / f"magma_v3_{name}.xml").write_bytes(raw)
        root = ET.fromstring(raw)
        warning = " ".join(
            "".join(w.itertext()).strip() for w in root.findall(".//warning")
        ).strip()
        lines = ["".join(line.itertext()) for line in root.findall(".//results/line")]
        output = "\n".join(lines)
        (RESULTS / f"magma_v3_{name}.txt").write_text(output + "\n", encoding="utf-8")
        error_markers = (
            "Runtime error", "User error", "Internal error", "Assertion failed",
            "Identifier '", "Syntax error",
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
            "output_file": f"magma_v3_{name}.txt",
            "xml_file": f"magma_v3_{name}.xml",
            "output_tail": output[-2000:],
        }
    except Exception:
        if raw:
            (RESULTS / f"magma_v3_{name}.xml").write_bytes(raw)
        return {
            "name": name,
            "status": "exception",
            "elapsed_seconds": time.time() - started,
            "traceback": traceback.format_exc(),
        }


def parse_number_list(text: str) -> list[decimal.Decimal]:
    left = text.find("[")
    right = text.rfind("]")
    if left < 0 or right < left:
        raise ValueError(f"missing list in {text[:120]!r}")
    values = [v.strip() for v in text[left + 1:right].split(",")]
    return [decimal.Decimal(v) for v in values if v]


def reconstruct_height_matrix(task_results: list[dict]) -> dict:
    decimal.getcontext().prec = 140
    rows: dict[int, list[decimal.Decimal]] = {}
    for result in task_results:
        if not result["name"].startswith("height_rows_80_") or result["status"] != "pass":
            continue
        text = (RESULTS / result["output_file"]).read_text(encoding="utf-8")
        for line in text.splitlines():
            match = re.match(r"^R30_HEIGHT_ROW\s+(\d+)\s+(\[.*\])$", line.strip())
            if match:
                rows[int(match.group(1))] = parse_number_list(match.group(2))
    if sorted(rows) != list(range(1, 30)):
        raise RuntimeError(f"height rows incomplete: {sorted(rows)}")
    if any(len(row) != 29 for row in rows.values()):
        raise RuntimeError("height row length is not 29")
    matrix = [rows[i] for i in range(1, 30)]
    max_symmetry_error = max(
        abs(matrix[i][j] - matrix[j][i]) for i in range(29) for j in range(29)
    )

    # High-precision Gaussian elimination on the printed 80-digit matrix.
    work = [row[:] for row in matrix]
    determinant = decimal.Decimal(1)
    sign = 1
    for col in range(29):
        pivot = max(range(col, 29), key=lambda r: abs(work[r][col]))
        if work[pivot][col] == 0:
            raise RuntimeError("printed height matrix is singular")
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            sign *= -1
        piv = work[col][col]
        determinant *= piv
        for row in range(col + 1, 29):
            factor = work[row][col] / piv
            if factor:
                for j in range(col + 1, 29):
                    work[row][j] -= factor * work[col][j]
            work[row][col] = decimal.Decimal(0)
    if sign < 0:
        determinant = -determinant

    summaries: dict[int, dict] = {}
    for precision in (80, 160):
        result = next(r for r in task_results if r["name"] == f"height_summary_{precision}")
        text = (RESULTS / result["output_file"]).read_text(encoding="utf-8")
        regulator_match = re.search(r"^R30_REGULATOR\s+(.+)$", text, re.MULTILINE)
        diagonal_match = re.search(r"^R30_HEIGHT_DIAGONAL\s+(\[.*\])$", text, re.MULTILINE)
        if not regulator_match or not diagonal_match:
            raise RuntimeError(f"cannot parse height summary {precision}")
        summaries[precision] = {
            "regulator": decimal.Decimal(regulator_match.group(1).strip()),
            "diagonal": parse_number_list(diagonal_match.group(1)),
        }
    max_diagonal_error_80 = max(
        abs(matrix[i][i] - summaries[80]["diagonal"][i]) for i in range(29)
    )
    relative_det_error = abs(determinant - summaries[80]["regulator"]) / abs(summaries[80]["regulator"])
    relative_summary_change = abs(summaries[160]["regulator"] - summaries[80]["regulator"]) / abs(summaries[160]["regulator"])

    certificate = {
        "precision_of_rows": 80,
        "matrix": [[str(v) for v in row] for row in matrix],
        "regulator_recomputed_from_printed_matrix": str(determinant),
        "regulator_reported_80": str(summaries[80]["regulator"]),
        "regulator_reported_160": str(summaries[160]["regulator"]),
        "maximum_symmetry_error": str(max_symmetry_error),
        "maximum_diagonal_error_against_summary_80": str(max_diagonal_error_80),
        "relative_determinant_error_reconstruction": str(relative_det_error),
        "relative_regulator_change_80_to_160": str(relative_summary_change),
        "checks": {
            "all_29_rows_present": True,
            "matrix_is_symmetric_to_printed_precision": max_symmetry_error < decimal.Decimal("1e-70"),
            "diagonal_matches_summary": max_diagonal_error_80 < decimal.Decimal("1e-70"),
            "recomputed_determinant_matches_summary": relative_det_error < decimal.Decimal("1e-65"),
            "regulator_stable_80_to_160": relative_summary_change < decimal.Decimal("1e-70"),
            "reported_regulators_nonzero": summaries[80]["regulator"] != 0 and summaries[160]["regulator"] != 0,
        },
    }
    certificate["pass"] = all(certificate["checks"].values())
    (RESULTS / "magma_v3_height_matrix.json").write_text(
        json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
    )
    with (RESULTS / "magma_v3_height_matrix.txt").open("w", encoding="utf-8") as fh:
        fh.write(f"# Magma V2.29, printed precision 80\n")
        fh.write(f"# regulator_80={summaries[80]['regulator']}\n")
        fh.write(f"# regulator_160={summaries[160]['regulator']}\n")
        for row in matrix:
            fh.write(" ".join(str(v) for v in row) + "\n")
    return certificate


def main() -> int:
    blocks = [(1, 5), (6, 10), (11, 15), (16, 20), (21, 25), (26, 29)]
    tasks: list[tuple[str, str, str]] = [
        ("model_torsion", task_model_torsion(), "R30_MODEL_TORSION_PASS"),
        ("local_mod2", task_local_mod2(), "R30_EXACT_LOCAL_INDEPENDENCE_PASS"),
        ("independence", task_independence(), "R30_DIRECT_INDEPENDENCE_PASS"),
        ("height_summary_80", task_height_summary(80), "R30_HEIGHT_SUMMARY_80_PASS"),
        ("height_summary_160", task_height_summary(160), "R30_HEIGHT_SUMMARY_160_PASS"),
    ]
    for first, last in blocks:
        name = f"height_rows_80_{first}_{last}"
        tasks.append((name, task_height_rows(80, first, last), f"R30_HEIGHT_ROWS_80_{first}_{last}_PASS"))

    results = []
    for name, code, sentinel in tasks:
        (RESULTS / f"magma_v3_{name}.m").write_text(code, encoding="utf-8")
        result = submit(name, code, sentinel)
        results.append(result)
        print(json.dumps(result, indent=2), flush=True)

    tasks_pass = all(r["status"] == "pass" for r in results)
    height_certificate: dict | None = None
    height_error: str | None = None
    if tasks_pass:
        try:
            height_certificate = reconstruct_height_matrix(results)
        except Exception:
            height_error = traceback.format_exc()
    core_pass = tasks_pass and bool(height_certificate and height_certificate.get("pass"))
    document = {
        "service": SERVER,
        "strict_parser": True,
        "calculator_limits_handled_by_row_blocks": True,
        "tasks": results,
        "height_certificate_file": "magma_v3_height_matrix.json" if height_certificate else None,
        "height_reconstruction_error": height_error,
        "core_pass": core_pass,
    }
    (RESULTS / "magma_free_result_v3.json").write_text(
        json.dumps(document, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"core_pass": core_pass, "height_error": height_error}, indent=2), flush=True)
    return 0 if core_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
