#!/usr/bin/env python3
"""Dipole probe DFT reflection for grading / SMA sweeps.

For each PointProbe (default component Ez for a z-directed dipole):
  - Incident window: around the first outgoing peak.
  - Reflection window: after the field drops below a fraction of the peak
    (quiet gap), through final_time — late return from outer truncation.

Reports 20*log10(|E_ref|/|E_inc|) at the incident spectral peak.

Example:
  python3 scripts/dipole_pml_grading_dft.py \\
    Exports/mpi-4/2D_Dipole Exports/mpi-4/2D_Dipole_PML_GO0 ...
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pml_dft_reflection import (  # noqa: E402
    C_SI,
    dft_spectrum,
    find_f_peak,
    parse_probe_dat,
    reflection_db,
    resolve_probe_path,
)

PROBE_LABELS = {0: "X", 1: "Y", 2: "XY"}


def auto_windows(
    t: np.ndarray,
    y: np.ndarray,
    *,
    quiet_frac: float,
    quiet_hold: float,
    ref_delay: float,
) -> tuple[tuple[float, float], tuple[float, float], dict]:
    """Return (inc_window, ref_window, info) in the same time unit as t.

    ref_delay: minimum code-time after t_peak before the reflection window may
    start (covers round-trip to the outer truncation).
    """
    abs_y = np.abs(y)
    i_peak = int(np.argmax(abs_y))
    peak = float(abs_y[i_peak])
    t_peak = float(t[i_peak])
    if peak <= 0.0:
        raise ValueError("zero peak in probe trace")

    thresh_rise = 0.05 * peak
    i0 = 0
    for i in range(i_peak, -1, -1):
        if abs_y[i] < thresh_rise:
            i0 = i
            break
    t_inc0 = float(t[i0])
    half = max(t_peak - t_inc0, 0.5)
    t_inc1 = min(float(t[-1]), t_peak + 1.25 * half)

    # Earliest allowed reflection start: after round-trip delay from peak.
    t_ref_earliest = t_peak + ref_delay

    quiet_level = quiet_frac * peak
    dt = float(np.median(np.diff(t)))
    need = max(1, int(round(quiet_hold / max(dt, 1e-12))))
    i_start = int(np.searchsorted(t, max(t_inc1, t_ref_earliest), side="left"))
    i_quiet = None
    run = 0
    for i in range(i_start, len(t)):
        if abs_y[i] < quiet_level:
            run += 1
            if run >= need:
                i_quiet = i - need + 1
                break
        else:
            run = 0

    if i_quiet is None:
        t_ref0 = t_ref_earliest
    else:
        t_ref0 = max(t_ref_earliest, float(t[i_quiet]))
    t_ref1 = float(t[-1])
    if t_ref1 - t_ref0 < 1.0:
        t_ref0 = max(t_inc1, t_ref1 - 2.0)

    info = {
        "t_peak": t_peak,
        "peak": peak,
        "t_inc": (t_inc0, t_inc1),
        "t_ref": (t_ref0, t_ref1),
        "quiet_level": quiet_level,
    }
    return (t_inc0, t_inc1), (t_ref0, t_ref1), info


def analyze_export(
    export_dir: Path,
    *,
    quiet_frac: float,
    quiet_hold: float,
    ref_delay: float,
    f_band: tuple[float, float],
    component: str,
) -> list[dict]:
    rows = []
    for pid in range(3):
        path = resolve_probe_path(export_dir, pid)
        t_si, fields, meta = parse_probe_dat(path)
        t = t_si * C_SI
        pos = meta.get("position") or (0.0, 0.0, 0.0)
        if component not in fields:
            raise KeyError(f"missing component {component}")
        y = fields[component]
        inc_w, ref_w, info = auto_windows(
            t, y, quiet_frac=quiet_frac, quiet_hold=quiet_hold, ref_delay=ref_delay
        )
        m_inc = (t >= inc_w[0]) & (t <= inc_w[1])
        m_ref = (t >= ref_w[0]) & (t <= ref_w[1])
        if not np.any(m_inc) or not np.any(m_ref):
            raise ValueError(f"{export_dir} probe {pid}: empty window")

        f_inc, a_inc = dft_spectrum(t[m_inc], y[m_inc], taper=True)
        f_ref, a_ref = dft_spectrum(t[m_ref], y[m_ref], taper=True)
        a_ref_i = np.interp(f_inc, f_ref, a_ref, left=0.0, right=0.0)
        r_db = reflection_db(a_ref_i, a_inc)
        fmin, fmax = f_band
        i_peak = find_f_peak(f_inc, a_inc, fmin, fmax)
        band = (f_inc >= fmin) & (f_inc <= fmax) & (f_inc > 0.0)
        i_worst = (
            int(np.flatnonzero(band)[int(np.argmax(r_db[band]))])
            if np.any(band)
            else i_peak
        )
        rows.append(
            {
                "case": export_dir.name,
                "probe": pid,
                "label": PROBE_LABELS.get(pid, str(pid)),
                "x": float(pos[0]),
                "y": float(pos[1]),
                "t_peak": info["t_peak"],
                "inc": inc_w,
                "ref": ref_w,
                "f_peak": float(f_inc[i_peak]),
                "R_dB": float(r_db[i_peak]),
                "R_worst": float(r_db[i_worst]),
                "max_inc": float(np.max(np.abs(y[m_inc]))),
                "max_ref": float(np.max(np.abs(y[m_ref]))),
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("exports", nargs="+", type=Path, help="Exports/.../<case> dirs")
    ap.add_argument("--quiet-frac", type=float, default=0.02)
    ap.add_argument("--quiet-hold", type=float, default=0.75, help="code-time")
    ap.add_argument(
        "--ref-delay",
        type=float,
        default=3.0,
        help="min code-time after peak before reflection window (default 3 ≈ "
        "round-trip from r=2.5 to outer wall at 4)",
    )
    ap.add_argument("--f-band", type=float, nargs=2, default=(0.05, 1.5))
    ap.add_argument(
        "--component",
        default="Ez",
        choices=("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"),
        help="Field component (default Ez for z-directed dipole)",
    )
    args = ap.parse_args()

    all_rows: list[dict] = []
    for d in args.exports:
        all_rows.extend(
            analyze_export(
                d,
                quiet_frac=args.quiet_frac,
                quiet_hold=args.quiet_hold,
                ref_delay=args.ref_delay,
                f_band=tuple(args.f_band),
                component=args.component,
            )
        )

    print(f"component: {args.component}")
    print(
        f"{'case':<22} {'probe':<4} {'f_peak':>8} {'R_dB':>8} {'R_worst':>8} "
        f"{'max|inc|':>10} {'max|ref|':>10}  inc_win  ref_win"
    )
    for r in all_rows:
        print(
            f"{r['case']:<22} {r['label']:<4} {r['f_peak']:8.4f} {r['R_dB']:8.2f} "
            f"{r['R_worst']:8.2f} {r['max_inc']:10.3e} {r['max_ref']:10.3e}  "
            f"[{r['inc'][0]:.2f},{r['inc'][1]:.2f}] [{r['ref'][0]:.2f},{r['ref'][1]:.2f}]"
        )

    print(f"\n--- summary @ X probe (R_dB at f_peak, {args.component}) ---")
    for r in all_rows:
        if r["label"] == "X":
            print(f"  {r['case']:<22}  {r['R_dB']:7.2f} dB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
