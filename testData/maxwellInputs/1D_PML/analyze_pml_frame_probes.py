#!/usr/bin/env python3
"""Correlate ParaView cycles with point probes (1D PML cases)."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

C_SI = 299792458.0
FINAL_CODE = 20.0
N_CYCLES = 199
T_NS = 66.713


def cycle_to_mult(cycle: int, steps: int = 20000, rk4: int = 4) -> int:
    step = round(cycle / N_CYCLES * steps)
    return step * rk4


def cycle_to_t_ns(cycle: int) -> float:
    return T_NS * cycle / N_CYCLES


def read_probe(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    with path.open() as f:
        for line in f:
            if line.startswith(("Point", "Spatial", "Time")) or not line.strip():
                continue
            p = line.split()
            if len(p) < 7:
                continue
            rows.append((float(p[0]), float(p[2]), float(p[6])))
    return rows


def sample_at_cycle(rows: list, cycle: int) -> tuple[float, float, float]:
    tp = (cycle / N_CYCLES * FINAL_CODE) / C_SI
    return min(rows, key=lambda r: abs(r[0] - tp))


def first_gate_break(rows: list, thr: float = 3.0) -> int | None:
    for i, r in enumerate(rows):
        if abs(r[1]) > thr or abs(r[2]) > thr:
            return i
    return None


def parse_mult_diag_log(path: Path) -> dict[int, dict]:
    pat = re.compile(
        r"\[PML Mult diag @ call (\d+)\] max\|psi\|=([\d.e+-]+) "
        r"outer=([\d.e+-]+) inner_pml=([\d.e+-]+) max\|dpsi/dt\|=([\d.e+-]+) "
        r"max\|out E\|=([\d.e+-]+) max\|out H\|=([\d.e+-]+) "
        r"max\|Ey\|=([\d.e+-]+) outer_Ey=([\d.e+-]+)"
    )
    out: dict[int, dict] = {}
    if not path.is_file():
        return out
    for line in path.read_text().splitlines():
        m = pat.search(line)
        if m:
            call = int(m.group(1))
            out[call] = {
                "max_psi": float(m.group(2)),
                "psi_outer": float(m.group(3)),
                "psi_inner": float(m.group(4)),
                "max_dpsi": float(m.group(5)),
                "max_out_e": float(m.group(6)),
                "max_out_h": float(m.group(7)),
                "max_ey": float(m.group(8)),
                "outer_ey": float(m.group(9)),
            }
    return out


def nearest_diag(diag: dict[int, dict], target_call: int) -> tuple[int, dict] | None:
    if not diag:
        return None
    call = min(diag.keys(), key=lambda c: abs(c - target_call))
    return call, diag[call]


def analyze_case(name: str, export_dir: Path, log_path: Path | None, cycles: list[int]) -> None:
    print(f"\n{'=' * 60}\n{name}\n{'=' * 60}")
    probes = {}
    for pid in range(3):
        p = export_dir / "PointProbes" / f"PointProbe{pid}.dat"
        probes[pid] = read_probe(p) if p.is_file() else []
        if probes[pid]:
            mx = max(max(abs(r[1]), abs(r[2])) for r in probes[pid])
            print(f"Probe {pid}: n={len(probes[pid])} max(|Ey|,|Hz|)={mx:.4e}")

    print(f"\n{'cyc':>4} {'t_ns':>7} {'mult':>7} | P0 Ey | P1 Ey | P2 Ey")
    for k in cycles:
        vals = [sample_at_cycle(probes[p], k) if probes[p] else (0, 0, 0) for p in range(3)]
        print(
            f"{k:4d} {cycle_to_t_ns(k):7.2f} {cycle_to_mult(k):7d} "
            f"{abs(vals[0][1]):9.4f} {abs(vals[1][1]):9.4f} {abs(vals[2][1]):9.4e}"
        )

    for pid, label in [(0, "1.99 vac"), (1, "2.01 PML")]:
        rows = probes.get(pid, [])
        if not rows:
            continue
        i = first_gate_break(rows)
        if i is not None:
            r = rows[i]
            tc = r[0] * C_SI
            print(
                f"Probe {pid} ({label}): first |E/H|>3 at row {i}, "
                f"t_ns={r[0]*1e9:.2f}, cycle~{round(tc/FINAL_CODE*N_CYCLES)}"
            )

    diag = parse_mult_diag_log(log_path) if log_path else {}
    if diag:
        print(f"\nPML Mult diag (from {log_path.name}):")
        print(f"{'cyc':>4} {'mult_tgt':>8} {'mult_log':>8} {'max|psi|':>12} {'outer':>10} {'outer_Ey':>10}")
        for k in cycles:
            tgt = cycle_to_mult(k)
            hit = nearest_diag(diag, tgt)
            if hit:
                c, d = hit
                print(
                    f"{k:4d} {tgt:8d} {c:8d} {d['max_psi']:12.4e} "
                    f"{d['psi_outer']:10.4e} {d['outer_ey']:10.4e}"
                )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument(
        "--cycles", type=int, nargs="+", default=[48, 61, 68, 75, 102, 144, 167, 199]
    )
    args = ap.parse_args()
    base: Path = args.base

    cases = [
        ("1D_PML_centered (α=0)", base / "1D_PML_centered", base / "1D_PML_centered" / "pml_diag.log"),
        ("1D_PML (α=1)", base / "1D_PML", base / "1D_PML" / "pml_diag.log"),
    ]
    for title, case_dir, log in cases:
        export_dir = case_dir / "Exports" / "single-core" / case_dir.name
        analyze_case(title, export_dir, log, args.cycles)
    return 0


if __name__ == "__main__":
    sys.exit(main())
