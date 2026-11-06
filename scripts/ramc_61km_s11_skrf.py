#!/usr/bin/env python3
"""Compare 1D RAM-C II 61 km stagnation-line S11 (90-110 GHz) from opensemba_dgtd against scikit-rf.

Setup (testData/testInputs/1D_Drude_RAMC_61km/1D_Drude_RAMC_61km.json):
  - Domain: [0, 1.0] m, with vacuum on [0, 0.988] m (0.5 mm elements, LGL order 4)
    and the 12 mm RAM-C II 61 km stagnation line on [0.988, 1.0] m (0.2 mm elements, LGL order 4).
  - Left boundary (x = 0.0 m): SMA; Right boundary (x = 1.0 m): PEC.
  - TFSF planewave entry at x = 0.25 m (+x propagation, modulated Gaussian centered at 100 GHz, f_1e = 20 GHz).
  - PointProbe 0 at x = 0.10 m (scattered-field region: pure reflected wave).
  - PointProbe 1 at x = 0.30 m (total-field region: incoming wave first, followed
    by reflection; cut at t_code = 1.00 m after the incoming wave has passed).
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import skrf as rf
from skrf.media import Freespace

C_SI = 299792458.0
REPO = Path(__file__).resolve().parents[1]
SOLVER = REPO / "build/gnu-release-mpi/bin/opensemba_dgtd"
CASE_JSON = REPO / "testData/testInputs/1D_Drude_RAMC_61km/1D_Drude_RAMC_61km.json"

X_REF = 0.10
X_INC = 0.30
X_INTERFACE = 0.988
T_CUT_CODE = 1.00  # light-metres (~3.34 ns): after incoming pulse at x=0.30 m, before reflection returns


def omega_p_ramc(x_local: np.ndarray | float) -> np.ndarray:
    x = np.asarray(x_local, dtype=float)
    wp_inf = 9.697700e+08
    x_sh = -7.309150e-03
    d_sh = 3.211841e-04
    a0 = 6.865065e+11
    a1 = -4.836413e+13
    a2 = -2.041384e+15
    a3 = 8.678069e+17
    d_bl = 4.296268e-04
    A_bl = 2.639789e+11
    step = 0.5 * (1.0 + np.tanh((x - x_sh) / d_sh))
    s = x + 0.0066
    env = a0 + s * (a1 + s * (a2 + s * a3)) - A_bl * np.exp(x / d_bl)
    return np.where(x > 0.0, 0.0, wp_inf + step * np.maximum(0.0, env))


def nu_c_gupta(x_local: np.ndarray | float) -> np.ndarray:
    x = np.asarray(x_local, dtype=float)
    nu_inf = 2.732100e+07
    x_sh = -7.482650e-03
    d_sh = 4.340546e-04
    a0 = 5.765874e+09
    a1 = 6.977597e+11
    a2 = -1.255220e+14
    a3 = 8.677992e+15
    d_bl = 2.616807e-04
    A_bl = 2.012032e+09
    step = 0.5 * (1.0 + np.tanh((x - x_sh) / d_sh))
    s = x + 0.0066
    env = a0 + s * (a1 + s * (a2 + s * a3)) + A_bl * np.exp(x / d_bl)
    return np.where(x > 0.0, 0.0, nu_inf + step * np.maximum(0.0, env))


def load_probe(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 7:
            continue
        try:
            rows.append([float(v) for v in parts])
        except ValueError:
            continue
    arr = np.asarray(rows)
    return arr[:, 0], arr[:, 2]  # t_si [s], Ey


def compute_skrf_s11(freqs_hz: np.ndarray, n_slices: int = 2400) -> np.ndarray:
    freq_obj = rf.Frequency.from_f(freqs_hz, unit="Hz")
    omega = 2.0 * np.pi * freqs_hz
    air = Freespace(frequency=freq_obj)

    edges = np.linspace(-0.012, 0.0, n_slices + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    dx = edges[1] - edges[0]

    net = air.thru()
    for xc in centers:
        wp = float(omega_p_ramc(xc))
        nu = float(nu_c_gupta(xc))
        eps_r = 1.0 - (wp**2) / (omega * (omega - 1j * nu))
        layer = Freespace(frequency=freq_obj, ep_r=eps_r, mu_r=1.0).line(d=dx, unit="m")
        net = net ** layer
    net = net ** air.short()
    return net.s[:, 0, 0]


def dft(time_s: np.ndarray, signal: np.ndarray, freqs_hz: np.ndarray) -> np.ndarray:
    dt = float(np.median(np.diff(time_s)))
    phase = np.exp(-2j * np.pi * freqs_hz[:, None] * time_s[None, :])
    return (phase * signal[None, :]).sum(axis=1) * dt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--np", type=int, default=8, help="Number of MPI ranks (default: 8)")
    parser.add_argument("--skip-run", action="store_true", help="Use existing probe outputs without re-running solver")
    args = parser.parse_args()

    run_mode = "single-core" if args.np == 1 else f"mpi-{args.np}"
    out_dir = REPO / "exports/SimulationData" / run_mode / "1D_Drude_RAMC_61km"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_run:
        log_path = out_dir / "solver.log"
        with log_path.open("w") as log_file:
            subprocess.run(
                ["mpirun", "-np", str(args.np), str(SOLVER), "-i", str(CASE_JSON), "--restart"],
                cwd=REPO,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,
            )

    t_ref, ey_ref = load_probe(out_dir / "PointProbes/PointProbe0.dat")
    t_inc, ey_inc_raw = load_probe(out_dir / "PointProbes/PointProbe1.dat")

    t_code = t_inc * C_SI
    ey_inc_cut = np.where(t_code <= T_CUT_CODE, ey_inc_raw, 0.0)

    # Frequency grid from 1 MHz (0.001 GHz) to 250 GHz
    freqs_hz = np.concatenate([
        np.geomspace(1.0e6, 1.0e9, 40, endpoint=False),
        np.linspace(1.0e9, 250.0e9, 997),
    ])
    F_inc = dft(t_inc, ey_inc_cut, freqs_hz)
    F_ref = dft(t_ref, ey_ref, freqs_hz)

    # De-embed phase from probe locations (x_inc=0.30 m, x_ref=0.10 m) to plasma interface x_int=0.988 m
    k0 = 2.0 * np.pi * freqs_hz / C_SI
    s11_dgtd = (F_ref / F_inc) * np.exp(1j * k0 * ((X_INTERFACE - X_REF) + (X_INTERFACE - X_INC)))
    s11_dgtd_db = 20.0 * np.log10(np.abs(s11_dgtd))
    s11_dgtd_deg = np.rad2deg(np.unwrap(np.angle(s11_dgtd)))

    s11_skrf = compute_skrf_s11(freqs_hz, n_slices=2400)
    s11_skrf_db = 20.0 * np.log10(np.abs(s11_skrf))
    s11_skrf_deg = np.rad2deg(np.unwrap(np.angle(s11_skrf)))

    csv_path = out_dir / "ramc_61km_1mhz_250ghz_s11.csv"
    with csv_path.open("w") as f:
        f.write("freq_hz,s11_dgtd_db,s11_skrf_db,s11_dgtd_deg_unwrapped,s11_skrf_deg_unwrapped\n")
        for freq, db_d, db_s, ph_d, ph_s in zip(
            freqs_hz, s11_dgtd_db, s11_skrf_db, s11_dgtd_deg, s11_skrf_deg
        ):
            f.write(f"{freq:.8e},{db_d:.6f},{db_s:.6f},{ph_d:.6f},{ph_s:.6f}\n")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))

    # Panel 1: |S11| in dB (1 MHz - 250 GHz)
    ax = axes[0]
    ax.plot(freqs_hz / 1e9, s11_skrf_db, color="black", linewidth=2.0, label="scikit-rf (2400 layers + PEC)")
    ax.plot(freqs_hz / 1e9, s11_dgtd_db, color="tab:red", linestyle="--", linewidth=1.8, label=f"opensemba_dgtd (order 4, {run_mode})")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel(r"$|S_{11}|$ [dB]")
    ax.set_title(r"$|S_{11}|$ Magnitude (1 MHz – 250 GHz)")
    ax.set_xlim(0.0, 250.0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5)

    # Panel 2: Unwrapped phase of S11 (1 MHz - 250 GHz)
    ax = axes[1]
    ax.plot(freqs_hz / 1e9, s11_skrf_deg, color="black", linewidth=2.0, label="scikit-rf (2400 layers + PEC)")
    ax.plot(freqs_hz / 1e9, s11_dgtd_deg, color="tab:red", linestyle="--", linewidth=1.8, label="opensemba_dgtd (de-embedded to 0.988 m)")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel(r"Unwrapped $\angle S_{11}$ [deg]")
    ax.set_title(r"$\angle S_{11}$ Phase at Plasma Interface ($x=0.988$ m)")
    ax.set_xlim(0.0, 250.0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5)

    fig.suptitle("RAM-C II 61 km Stagnation Line (1 MHz – 250 GHz, Order 4) — DGTD vs scikit-rf", fontsize=12)
    fig.tight_layout()
    png_path = out_dir / "ramc_61km_1mhz_250ghz_s11.png"
    fig.savefig(png_path, dpi=150)

    err_db = np.abs(s11_dgtd_db - s11_skrf_db)
    err_deg = np.abs((s11_dgtd_deg - s11_skrf_deg + 180.0) % 360.0 - 180.0)

    print(f"Wrote {csv_path}")
    print(f"Wrote {png_path}")
    print(f"Band [1 MHz, 250 GHz]: mean |ΔS11| = {np.mean(err_db):.4f} dB, max |ΔS11| = {np.max(err_db):.4f} dB, mean |Δ∠S11| = {np.mean(err_deg):.3f} deg, max |Δ∠S11| = {np.max(err_deg):.3f} deg")


if __name__ == "__main__":
    main()
