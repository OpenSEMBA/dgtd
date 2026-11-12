#!/usr/bin/env python3
"""Compare 1D RAM-C II 50 km stagnation-line S11 against scikit-rf.

Setup (testData/testInputs/1D_Drude_RAMC_50km/1D_Drude_RAMC_50km.json):
  - Domain [0, 1] m. Vacuum on [0, 0.976] m (1 mm elements).
    The 24 mm RAM-C II 50 km stagnation line sits on [0.976, 1] m
    (0.0125 mm elements, LGL order 2).
  - Left boundary x = 0: SMA. Right boundary x = 1: PEC.
  - TFSF planewave at x = 0.25 m, +x, Gaussian with f_1e = 5 GHz.
  - PointProbe 0 at x = 0.10 m (scattered field, reflected wave only).
  - PointProbe 1 at x = 0.30 m (total field). Incoming samples are kept
    until t_code = 0.50 m, before the reflection returns.
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
CASE_JSON = REPO / "testData/testInputs/1D_Drude_RAMC_50km/1D_Drude_RAMC_50km.json"

X_REF = 0.10
X_INC = 0.30
X_INTERFACE = 0.976
T_CUT_CODE = 0.50


def omega_p_ramc(x_local: np.ndarray | float) -> np.ndarray:
    x = np.asarray(x_local, dtype=float)
    wp_inf = 1.891488933e+05
    x_sh = -1.404091339e-02
    d_sh = 3.346945412e-03
    a0 = -1.442830910e+08
    a1 = -8.139327178e+09
    a2 = 2.023092321e+12
    a3 = 1.607015069e+14
    d_bl = 3.531895113e-03
    A_bl = 5.860432868e+09
    step = 0.5 * (1.0 + np.tanh((x - x_sh) / d_sh))
    s = x + 1.018700058e-02
    env = a0 + s * (a1 + s * (a2 + s * a3)) + A_bl * np.exp(x / d_bl)
    return np.where(x > 0.0, 0.0, wp_inf + step * np.maximum(0.0, env))


def nu_c_gupta(x_local: np.ndarray | float) -> np.ndarray:
    x = np.asarray(x_local, dtype=float)
    nu_inf = 1.007955905e+08
    x_sh = -1.878338741e-02
    d_sh = 1.433905287e-03
    a0 = 5.432719481e+09
    a1 = 3.957922427e+10
    a2 = -4.318630790e+12
    a3 = 3.187292185e+14
    d_bl = 4.727815008e-04
    A_bl = 2.680914883e+09
    step = 0.5 * (1.0 + np.tanh((x - x_sh) / d_sh))
    s = x - 6.151333725e-04
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
    return arr[:, 0], arr[:, 2]


def compute_skrf_s11(freqs_hz: np.ndarray, n_slices: int = 2400) -> np.ndarray:
    freq_obj = rf.Frequency.from_f(freqs_hz, unit="Hz")
    omega = 2.0 * np.pi * freqs_hz
    air = Freespace(frequency=freq_obj)

    edges = np.linspace(-0.024, 0.0, n_slices + 1)
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
    parser.add_argument("--np", type=int, default=1, help="MPI ranks (default: 1)")
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()

    run_mode = "single-core" if args.np == 1 else f"mpi-{args.np}"
    out_dir = REPO / "exports/SimulationData" / run_mode / "1D_Drude_RAMC_50km"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_run:
        log_path = out_dir / "solver.log"
        cmd = [str(SOLVER), "-i", str(CASE_JSON), "--restart"]
        if args.np > 1:
            cmd = ["mpirun", "-np", str(args.np), *cmd]
        with log_path.open("w") as log_file:
            subprocess.run(cmd, cwd=REPO, stdout=log_file, stderr=subprocess.STDOUT, check=True)

    t_ref, ey_ref = load_probe(out_dir / "PointProbes/PointProbe0.dat")
    t_inc, ey_inc_raw = load_probe(out_dir / "PointProbes/PointProbe1.dat")

    t_code = t_inc * C_SI
    ey_inc_cut = np.where(t_code <= T_CUT_CODE, ey_inc_raw, 0.0)

    freqs_hz = np.concatenate([
        np.geomspace(1.0e6, 1.0e8, 60, endpoint=False),
        np.linspace(1.0e8, 5.0e9, 490),
    ])
    F_inc = dft(t_inc, ey_inc_cut, freqs_hz)
    F_ref = dft(t_ref, ey_ref, freqs_hz)

    k0 = 2.0 * np.pi * freqs_hz / C_SI
    s11_dgtd = (F_ref / F_inc) * np.exp(1j * k0 * ((X_INTERFACE - X_REF) + (X_INTERFACE - X_INC)))
    s11_dgtd_db = 20.0 * np.log10(np.abs(s11_dgtd))
    s11_dgtd_deg = np.rad2deg(np.unwrap(np.angle(s11_dgtd)))

    s11_skrf = compute_skrf_s11(freqs_hz, n_slices=2400)
    s11_skrf_db = 20.0 * np.log10(np.abs(s11_skrf))
    s11_skrf_deg = np.rad2deg(np.unwrap(np.angle(s11_skrf)))

    csv_path = out_dir / "ramc_50km_1mhz_5ghz_s11.csv"
    with csv_path.open("w") as f:
        f.write("freq_hz,s11_dgtd_db,s11_skrf_db,s11_dgtd_deg_unwrapped,s11_skrf_deg_unwrapped\n")
        for freq, db_d, db_s, ph_d, ph_s in zip(
            freqs_hz, s11_dgtd_db, s11_skrf_db, s11_dgtd_deg, s11_skrf_deg
        ):
            f.write(f"{freq:.8e},{db_d:.6f},{db_s:.6f},{ph_d:.6f},{ph_s:.6f}\n")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    ax = axes[0]
    ax.plot(freqs_hz / 1e9, s11_skrf_db, color="black", linewidth=2.0, label="scikit-rf (2400 layers + PEC)")
    ax.plot(freqs_hz / 1e9, s11_dgtd_db, color="tab:red", linestyle="--", linewidth=1.8,
            label=f"opensemba_dgtd (LGL order 2, {run_mode})")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel(r"$|S_{11}|$ [dB]")
    ax.set_title(r"$|S_{11}|$ (1 MHz – 5 GHz)")
    ax.set_xlim(0.0, 5.0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5)

    ax = axes[1]
    ax.plot(freqs_hz / 1e9, s11_skrf_deg, color="black", linewidth=2.0, label="scikit-rf")
    ax.plot(freqs_hz / 1e9, s11_dgtd_deg, color="tab:red", linestyle="--", linewidth=1.8,
            label="opensemba_dgtd (de-embedded to 0.976 m)")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel(r"Unwrapped $\angle S_{11}$ [deg]")
    ax.set_title(r"$\angle S_{11}$ at the plasma face ($x=0.976$ m)")
    ax.set_xlim(0.0, 5.0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5)

    fig.suptitle("RAM-C II 50 km stagnation line (1 MHz – 5 GHz, LGL order 2)", fontsize=12)
    fig.tight_layout()
    png_path = out_dir / "ramc_50km_1mhz_5ghz_s11.png"
    fig.savefig(png_path, dpi=150)

    err_db = np.abs(s11_dgtd_db - s11_skrf_db)
    err_deg = np.abs((s11_dgtd_deg - s11_skrf_deg + 180.0) % 360.0 - 180.0)
    print(f"Wrote {csv_path}")
    print(f"Wrote {png_path}")
    print(
        "Band [1 MHz, 5 GHz]: "
        f"mean |ΔS11| = {np.mean(err_db):.4f} dB, max |ΔS11| = {np.max(err_db):.4f} dB, "
        f"mean |Δ∠S11| = {np.mean(err_deg):.3f} deg, max |Δ∠S11| = {np.max(err_deg):.3f} deg"
    )


if __name__ == "__main__":
    main()
