#!/usr/bin/env python3
"""Analyze probe matrix results for cross-configuration agreement."""
from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

SAME_CASE_RTOL = 0.05
PAIR_RTOL = 0.15
FIELD_MIN = 1e-12

PAIRS = [
    ("2D_RCS_Circle_G1", "2D_RCS_Circle_G1_hesthaven"),
    ("2D_RCS_Circle_G2", "2D_RCS_Circle_G2_hesthaven"),
    ("2D_RCS_PEC_Circle_G1_1m", "2D_RCS_PEC_Circle_G1_1m_Hesthaven"),
    ("3D_Nasa_Almond_G1", "3D_Nasa_Almond_G1_hesthaven"),
]


def mesh_name(case: str, root: Path) -> str:
    import json

    json_path = root / "testData/maxwellInputs" / case / f"{case}.json"
    data = json.loads(json_path.read_text())
    name = data["model"]["filename"]
    return name[:-4] if name.endswith(".msh") else name


def load_rows(path: Path):
    rows = []
    with path.open() as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            case, preset, np_s, status, elapsed, probe_path, max_e, final_ey = parts[:8]
            rows.append(
                {
                    "case": case,
                    "preset": preset,
                    "np": int(np_s),
                    "status": status,
                    "max_e": float(max_e),
                    "final_ey": float(final_ey),
                    "probe_path": probe_path,
                }
            )
    return rows


def rel_diff(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), 1e-30)
    return abs(a - b) / denom


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    results_path = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "build/probe_matrix_logs/results.tsv"
    rows = load_rows(results_path)

    print(f"\n=== Probe matrix analysis ({results_path}) ===\n")

    by_mesh: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row["status"] != "PASS":
            continue
        m = mesh_name(row["case"], root)
        by_mesh[m].append(row)

    exit_code = 0

    print("Per-case configuration agreement (max |E|):")
    for mesh, group in sorted(by_mesh.items()):
        if len(group) < 2:
            print(f"  {mesh}: only {len(group)} successful run(s), skip compare")
            continue
        ref = group[0]["max_e"]
        worst = 0.0
        worst_tag = ""
        for g in group[1:]:
            rd = rel_diff(g["max_e"], ref)
            if rd > worst:
                worst = rd
                worst_tag = f"{g['preset']} np{g['np']}"
        ok = worst <= SAME_CASE_RTOL
        mark = "OK" if ok else "FAIL"
        print(f"  [{mark}] {mesh}: worst rel diff {worst:.3e} vs ref ({worst_tag})")
        if not ok:
            exit_code = 1

    print("\nGlobal vs Hesthaven pairs (max |E|):")
    mesh_to_case = {mesh_name(c, root): c for c in {r["case"] for r in rows}}
    for global_case, hest_case in PAIRS:
        g_mesh = mesh_name(global_case, root)
        h_mesh = mesh_name(hest_case, root)
        g_vals = [r["max_e"] for r in by_mesh.get(g_mesh, [])]
        h_vals = [r["max_e"] for r in by_mesh.get(h_mesh, [])]
        if not g_vals or not h_vals:
            print(f"  [SKIP] {global_case} vs {hest_case}: missing successful runs")
            continue
        g_ref = sum(g_vals) / len(g_vals)
        h_ref = sum(h_vals) / len(h_vals)
        rd = rel_diff(g_ref, h_ref)
        ok = rd <= PAIR_RTOL
        mark = "OK" if ok else "FAIL"
        print(f"  [{mark}] {global_case} vs {hest_case}: rel diff {rd:.3e}")
        if not ok:
            exit_code = 1

    print("\nFailed / skipped runs:")
    any_bad = False
    for row in rows:
        if row["status"] != "PASS":
            any_bad = True
            print(f"  {row['case']} {row['preset']} np{row['np']}: {row['status']}")
    if not any_bad:
        print("  (none)")

    print("\nField detection (max |E| > {:.0e}):".format(FIELD_MIN))
    for row in rows:
        if row["status"] == "PASS" and row["max_e"] <= FIELD_MIN:
            print(f"  FAIL {row['case']} {row['preset']} np{row['np']}: max|E|={row['max_e']:.3e}")
            exit_code = 1

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
