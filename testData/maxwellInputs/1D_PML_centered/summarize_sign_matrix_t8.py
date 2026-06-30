#!/usr/bin/env python3
"""Summarize run_sign_matrix_t8.sh logs and per-mode probe exports."""
from __future__ import annotations

import re
import sys
from pathlib import Path

C_SI = 299792458.0
FINAL_CODE = 8.05
NAMES = {
    0: "default",
    1: "S1_flip_corrections",
    2: "S2_face_same_as_vol",
    3: "S3_face_cross_column",
    4: "S4_negate_driver_w",
    5: "S5_negate_pml_mult",
    6: "S6_flip_psi_mass",
}

DIAG = re.compile(
    r"\[PML Mult diag @ call (\d+)\] max\|psi\|=([\d.e+-]+) "
    r"outer=([\d.e+-]+) inner_pml=([\d.e+-]+) max\|dpsi/dt\|=([\d.e+-]+) "
    r"max\|out E\|=([\d.e+-]+) max\|out H\|=([\d.e+-]+) "
    r"max\|Ey\|=([\d.e+-]+) outer_Ey=([\d.e+-]+)"
)


def last_diag(log: Path) -> dict | None:
    last = None
    for line in log.read_text().splitlines():
        m = DIAG.search(line)
        if m:
            last = {
                "call": int(m.group(1)),
                "max_psi": float(m.group(2)),
                "outer_Ey": float(m.group(9)),
            }
    return last


def read_probe(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    with path.open() as f:
        for line in f:
            if line.startswith("Time (s)"):
                continue
            p = line.split()
            if len(p) < 7:
                continue
            try:
                rows.append((float(p[0]), abs(float(p[2])), abs(float(p[6]))))
            except ValueError:
                continue
    return rows


def probe_stats(export: Path) -> tuple[float, float, float, str]:
    base = export / "PointProbes"
    rows0 = read_probe(base / "PointProbe0.dat")
    rows1 = read_probe(base / "PointProbe1.dat")
    if not rows0 or not rows1:
        return float("nan"), float("nan"), float("nan"), "?"
    t_end = FINAL_CODE / C_SI
    e0 = min(rows0, key=lambda r: abs(r[0] - t_end))[1]
    e1 = min(rows1, key=lambda r: abs(r[0] - t_end))[1]
    max_ey = max(max(r[1] for r in rows0), max(r[1] for r in rows1))
    max_hz = max(max(r[2] for r in rows0), max(r[2] for r in rows1))
    gate = "PASS" if max_ey <= 3.0 and max_hz <= 3.0 else "FAIL"
    return e0, e1, max_ey, gate


def main() -> int:
    case = Path(__file__).resolve().parent
    out = case / "sign_matrix_t8"

    print(
        f"{'mode':<6} {'name':<24} {'max|psi|':>12} {'outer_Ey':>10} "
        f"{'P0|Ey|':>10} {'P1|Ey|':>10} {'max|Ey|':>10} {'gate':>5}"
    )
    print("-" * 100)

    best_psi = None
    best_mode = None
    rows_out: list[str] = []

    for mode in range(7):
        log = out / f"S{mode}.log"
        exp = out / f"exports_S{mode}"
        d = last_diag(log) if log.is_file() else None
        if d is None:
            print(f"{mode:<6} {NAMES[mode]:<24}  (no log)")
            continue
        p0, p1, max_ey, gate = probe_stats(exp) if exp.is_dir() else (
            float("nan"),
            float("nan"),
            float("nan"),
            "?",
        )
        line = (
            f"{mode:<6} {NAMES[mode]:<24} {d['max_psi']:12.4g} {d['outer_Ey']:10.4g} "
            f"{p0:10.4g} {p1:10.4g} {max_ey:10.4g} {gate:>5}"
        )
        print(line)
        rows_out.append(line)
        if best_psi is None or d["max_psi"] < best_psi:
            best_psi = d["max_psi"]
            best_mode = mode

    print()
    if best_mode is not None:
        print(
            f"Lowest final max|psi|: S{best_mode} ({NAMES[best_mode]}) = {best_psi:.4g} "
            f"(baseline default S0 for comparison)"
        )
    ref = last_diag(out / "S0.log")
    if ref and best_psi is not None and best_mode != 0:
        ratio = best_psi / ref["max_psi"] if ref["max_psi"] else float("inf")
        print(f"Best vs default max|psi| ratio: {ratio:.4g}")

    summary = out / "summary.txt"
    summary.write_text("\n".join(rows_out) + "\n")
    print(f"Wrote {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
