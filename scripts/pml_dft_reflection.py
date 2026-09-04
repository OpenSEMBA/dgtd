#!/usr/bin/env python3
"""Windowed DFT reflection coefficient from dgtd PointProbe exports.

Acceptance (docs/pml/05-verification.md): 20*log10(|E_ref|/|E_inc|) <= -40 dB
in the frequency domain, using a vacuum probe with time-separated incident and
reflected lobes.

Probe .dat time column is SI seconds (t_code / c_SI). This script converts to
normalized code time with c_SI = 299792458 unless --time-unit si is set.

Example:
  python3 scripts/pml_dft_reflection.py \\
    Exports/single-core/1D_PML_DFT \\
    --probe 0 --component Ey \\
    --inc-window 3.0 7.0 --ref-window 8.5 14.0
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

C_SI = 299792458.0
COLUMNS = ("time", "Ex", "Ey", "Ez", "Hx", "Hy", "Hz")
DEFAULT_GATE_DB = -40.0


def parse_probe_dat(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], dict]:
    """Return (t_si, fields, meta) from a PointProbeN.dat file."""
    meta: dict = {"path": str(path), "position": None, "probe_id": None}
    times: list[float] = []
    fields: dict[str, list[float]] = {c: [] for c in COLUMNS[1:]}
    in_data = False

    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("PointProbe ID"):
                try:
                    meta["probe_id"] = int(line.split()[-1])
                except ValueError:
                    pass
                continue
            if line.startswith("Spatial Position"):
                continue
            # Position line: three floats after the Spatial Position header
            if meta["position"] is None and not in_data:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        meta["position"] = tuple(float(x) for x in parts[:3])
                        continue
                    except ValueError:
                        pass
            if line.startswith("Time ("):
                in_data = True
                continue
            if not in_data:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            times.append(float(parts[0]))
            for i, name in enumerate(COLUMNS[1:]):
                fields[name].append(float(parts[i + 1]))

    if not times:
        raise ValueError(f"no data rows in {path}")

    t = np.asarray(times, dtype=float)
    arr = {k: np.asarray(v, dtype=float) for k, v in fields.items()}
    return t, arr, meta


def resolve_probe_path(export_dir: Path, probe_id: int) -> Path:
    base = export_dir
    if (base / "PointProbes").is_dir():
        base = base / "PointProbes"
    path = base / f"PointProbe{probe_id}.dat"
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def extract_window(
    t: np.ndarray, y: np.ndarray, t0: float, t1: float
) -> tuple[np.ndarray, np.ndarray]:
    mask = (t >= t0) & (t <= t1)
    if not np.any(mask):
        raise ValueError(f"empty window [{t0}, {t1}] (t span [{t[0]}, {t[-1]}])")
    return t[mask], y[mask]


def dft_spectrum(
    t: np.ndarray, y: np.ndarray, *, taper: bool
) -> tuple[np.ndarray, np.ndarray]:
    """Return (freq_code, |Y|) from a real signal sampled on possibly uneven t.

    Resamples onto a uniform grid spanning the window (linear interp) then rFFT.
    Frequencies are in normalized code units (1/t_code), same as dgtd.
    """
    if t.size < 4:
        raise ValueError(f"need at least 4 samples in window, got {t.size}")

    dt = float(np.median(np.diff(t)))
    if dt <= 0.0:
        raise ValueError("non-positive median dt in window")

    t_u = np.arange(t[0], t[-1] + 0.5 * dt, dt)
    y_u = np.interp(t_u, t, y)
    if taper:
        y_u = y_u * np.hanning(y_u.size)

    y_hat = np.fft.rfft(y_u)
    freq = np.fft.rfftfreq(y_u.size, d=dt)
    return freq, np.abs(y_hat)


def reflection_db(e_ref: np.ndarray, e_inc: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    ratio = e_ref / np.maximum(e_inc, eps)
    return 20.0 * np.log10(np.maximum(ratio, eps))


def find_f_peak(freq: np.ndarray, amp: np.ndarray, fmin: float, fmax: float) -> int:
    mask = (freq >= fmin) & (freq <= fmax) & (freq > 0.0)
    if not np.any(mask):
        # fall back to global positive-frequency peak
        mask = freq > 0.0
    idx_local = int(np.argmax(amp[mask]))
    return int(np.flatnonzero(mask)[idx_local])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "export_dir",
        type=Path,
        help="Exports/<mode>/<case> or .../PointProbes",
    )
    parser.add_argument("--probe", type=int, default=0, help="PointProbe ID (default 0)")
    parser.add_argument(
        "--component",
        default="Ey",
        choices=COLUMNS[1:],
        help="Field component (default Ey)",
    )
    parser.add_argument(
        "--inc-window",
        type=float,
        nargs=2,
        metavar=("T0", "T1"),
        default=(3.0, 7.2),
        help="Incident lobe window in code time (default 3 7.2; tuned for 1D_PML_DFT P0)",
    )
    parser.add_argument(
        "--ref-window",
        type=float,
        nargs=2,
        metavar=("T0", "T1"),
        default=(8.5, 11.5),
        help="Reflected lobe window in code time (default 8.5 11.5; tuned for 1D_PML_DFT P0)",
    )
    parser.add_argument(
        "--f-band",
        type=float,
        nargs=2,
        metavar=("FMIN", "FMAX"),
        default=(0.05, 1.5),
        help="Frequency band (code units) for peak-incident report (default 0.05 1.5)",
    )
    parser.add_argument(
        "--gate-db",
        type=float,
        default=DEFAULT_GATE_DB,
        help=f"Pass if R_dB <= this at f_peak (default {DEFAULT_GATE_DB})",
    )
    parser.add_argument(
        "--no-taper",
        action="store_true",
        help="Disable Hann taper before DFT",
    )
    parser.add_argument(
        "--time-unit",
        choices=("code", "si"),
        default="code",
        help="Unit of --inc-window/--ref-window (default code)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional path to write freq, |Einc|, |Eref|, R_dB",
    )
    args = parser.parse_args()

    path = resolve_probe_path(args.export_dir, args.probe)
    t_si, fields, meta = parse_probe_dat(path)
    t_code = t_si * C_SI

    if args.time_unit == "si":
        t_work = t_si
        unit = "SI s"
    else:
        t_work = t_code
        unit = "code"

    y = fields[args.component]
    t_inc, y_inc = extract_window(t_work, y, args.inc_window[0], args.inc_window[1])
    t_ref, y_ref = extract_window(t_work, y, args.ref_window[0], args.ref_window[1])

    # DFTs on code-time abscissa so frequencies are in code units either way.
    if args.time_unit == "si":
        t_inc_c = t_inc * C_SI
        t_ref_c = t_ref * C_SI
    else:
        t_inc_c, t_ref_c = t_inc, t_ref

    taper = not args.no_taper
    f_inc, a_inc = dft_spectrum(t_inc_c, y_inc, taper=taper)
    f_ref, a_ref = dft_spectrum(t_ref_c, y_ref, taper=taper)

    # Interpolate reflected spectrum onto incident frequency grid.
    a_ref_i = np.interp(f_inc, f_ref, a_ref, left=0.0, right=0.0)
    r_db = reflection_db(a_ref_i, a_inc)

    fmin, fmax = args.f_band
    i_peak = find_f_peak(f_inc, a_inc, fmin, fmax)
    f_peak = float(f_inc[i_peak])
    r_at_peak = float(r_db[i_peak])

    band = (f_inc >= fmin) & (f_inc <= fmax) & (f_inc > 0.0)
    if np.any(band):
        i_worst = int(np.flatnonzero(band)[int(np.argmax(r_db[band]))])
        f_worst = float(f_inc[i_worst])
        r_worst = float(r_db[i_worst])
    else:
        f_worst, r_worst = f_peak, r_at_peak

    passed = r_at_peak <= args.gate_db
    status = "PASS" if passed else "FAIL"

    pos = meta.get("position")
    pos_str = (
        f"({pos[0]:.4g}, {pos[1]:.4g}, {pos[2]:.4g})" if pos is not None else "?"
    )

    print(f"probe file: {path}")
    print(f"probe id:   {meta.get('probe_id')}  position: {pos_str}")
    print(f"component:  {args.component}")
    print(f"time unit:  {unit}  (export t_si -> t_code = t_si * c_SI)")
    print(
        f"inc window: [{args.inc_window[0]}, {args.inc_window[1]}] {unit}  "
        f"n={t_inc.size}  max|y|={np.max(np.abs(y_inc)):.6e}"
    )
    print(
        f"ref window: [{args.ref_window[0]}, {args.ref_window[1]}] {unit}  "
        f"n={t_ref.size}  max|y|={np.max(np.abs(y_ref)):.6e}"
    )
    print(f"taper:      {'hann' if taper else 'none'}")
    print(f"f_peak (max |E_inc| in [{fmin}, {fmax}]): {f_peak:.6g}")
    print(f"|E_inc|(f_peak)={a_inc[i_peak]:.6e}  |E_ref|(f_peak)={a_ref_i[i_peak]:.6e}")
    print(f"R_dB(f_peak) = {r_at_peak:.2f} dB")
    print(f"worst R_dB in band = {r_worst:.2f} dB at f={f_worst:.6g}")
    print(f"gate: R_dB(f_peak) <= {args.gate_db} dB  =>  {status}")

    if args.csv is not None:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w") as ofs:
            ofs.write("freq_code,abs_Einc,abs_Eref,R_dB\n")
            for f, ai, ar, rd in zip(f_inc, a_inc, a_ref_i, r_db):
                ofs.write(f"{f:.12e},{ai:.12e},{ar:.12e},{rd:.12e}\n")
        print(f"wrote {args.csv}")

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
