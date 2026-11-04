#!/usr/bin/env python3
"""Generate scikit-rf reference S-parameter CSVs for the 1D SGBC paper benchmarks.

Benchmarks from Section IV of docs/main.pdf:
  - IV-A: 1D_SGBC_Multilayer (PEC-backed 3-layer RAM coating, 4.5 to 7.4 GHz)
  - IV-B: 1D_SGBC_Resonant   (Passthrough 3-layer symmetric slab, 5.0 to 7.0 GHz)

The frequency grids match the DFT bins of each case's PointProbe time series:
  - 1D_SGBC_Multilayer: N = 12001 samples, dt_probe = 0.001 / c_0 (t_final = 12 m)
  - 1D_SGBC_Resonant:   N = 10002 samples, dt_probe = 0.0012 / c_0 (t_final = 12 m)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import skrf as rf
from skrf.media import Freespace

C_SI = 299792458.0
EPS0_SI = 8.8541878128e-12
REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "testData" / "sgbcReferences"


def probe_dft_frequencies(case_json_path: Path, fmin_hz: float, fmax_hz: float) -> np.ndarray:
    data = json.loads(case_json_path.read_text())
    dt_code = float(data["solver_options"]["time_step"])
    t_final = float(data["solver_options"]["final_time"])
    saves = int(data["probes"]["point"][0].get("saves", 0))

    total_steps = int(math.ceil(t_final / dt_code))
    vis_steps = max(1, total_steps // saves) if saves > 0 else 1

    # Simulate the exact cycle sampling in ProbesManager::updateProbe:
    # sampled when cycle % vis_steps == 0 or abs(time - final_time) < 1e-8
    n_samples = 0
    t = 0.0
    cycle = 0
    while True:
        if abs(t - t_final) < 1e-8 or (cycle % vis_steps == 0):
            n_samples += 1
        if t > t_final - 1e-8 * dt_code:
            break
        step_dt = min(dt_code, t_final - t)
        t += step_dt
        cycle += 1

    dt_probe_si = (vis_steps * dt_code) / C_SI
    freqs = np.fft.rfftfreq(n_samples, d=dt_probe_si)
    mask = (freqs >= fmin_hz) & (freqs <= fmax_hz)
    return freqs[mask]


def build_multilayer_network(freq_hz: np.ndarray, layers: list[dict], pec_backed: bool) -> rf.Network:
    freq_obj = rf.Frequency.from_f(freq_hz, unit="Hz")
    air = Freespace(frequency=freq_obj)
    omega = 2.0 * np.pi * freq_obj.f

    net = air.thru()
    for layer in layers:
        eps_r = float(layer.get("relative_permittivity", 1.0))
        mu_r = float(layer.get("relative_permeability", 1.0))
        sigma = float(layer.get("bulk_conductivity", 0.0))
        width = float(layer["material_width"])
        eps_complex = eps_r - 1j * sigma / (omega * EPS0_SI)
        medium = Freespace(frequency=freq_obj, ep_r=eps_complex, mu_r=mu_r)
        net = net ** medium.line(d=width, unit="m")

    if pec_backed:
        net = net ** air.short()
    else:
        net = net ** air.thru()
    return net


def generate_pec_ram_reference() -> Path:
    case_path = REPO / "testData" / "maxwellInputs" / "1D_SGBC_Multilayer" / "1D_SGBC_Multilayer.json"
    data = json.loads(case_path.read_text())
    sgbc = next(b for b in data["model"]["boundaries"] if b["type"] == "SGBC")
    freqs = probe_dft_frequencies(case_path, 4.5e9, 7.4e9)
    net = build_multilayer_network(freqs, sgbc["layers"], pec_backed=True)
    s11_db = 20.0 * np.log10(np.abs(net.s[:, 0, 0]))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = OUT_DIR / "1D_SGBC_Multilayer_skrf.csv"
    with out_csv.open("w") as f:
        f.write("freq_hz,s11_db\n")
        for freq, s11 in zip(freqs, s11_db):
            f.write(f"{freq:.12e},{s11:.12e}\n")
    print(f"Wrote {out_csv} ({len(freqs)} frequencies, [{freqs[0]/1e9:.4f}, {freqs[-1]/1e9:.4f}] GHz)")
    return out_csv


def generate_passthrough_reference() -> Path:
    case_path = REPO / "testData" / "maxwellInputs" / "1D_SGBC_Resonant" / "1D_SGBC_Resonant.json"
    data = json.loads(case_path.read_text())
    sgbc = next(b for b in data["model"]["boundaries"] if b["type"] == "SGBC")
    freqs = probe_dft_frequencies(case_path, 5.0e9, 7.0e9)
    net = build_multilayer_network(freqs, sgbc["layers"], pec_backed=False)
    s11_db = 20.0 * np.log10(np.abs(net.s[:, 0, 0]))
    s21_db = 20.0 * np.log10(np.abs(net.s[:, 1, 0]))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = OUT_DIR / "1D_SGBC_Resonant_skrf.csv"
    with out_csv.open("w") as f:
        f.write("freq_hz,s11_db,s21_db\n")
        for freq, s11, s21 in zip(freqs, s11_db, s21_db):
            f.write(f"{freq:.12e},{s11:.12e},{s21:.12e}\n")
    print(f"Wrote {out_csv} ({len(freqs)} frequencies, [{freqs[0]/1e9:.4f}, {freqs[-1]/1e9:.4f}] GHz)")
    return out_csv


def generate_through_fast_references() -> tuple[Path, Path]:
    case_path = REPO / "testData" / "maxwellInputs" / "1D_SGBC_Through" / "1D_SGBC_Through.json"
    data = json.loads(case_path.read_text())
    sgbc = next(b for b in data["model"]["boundaries"] if b["type"] == "SGBC")
    freqs = probe_dft_frequencies(case_path, 0.2e9, 1.5e9)

    # 1) Passthrough single lossy conductive sheet (eps_r=1, mu_r=1, sigma=100 S/m, d=5 mm)
    net_through = build_multilayer_network(freqs, [sgbc["material"]], pec_backed=False)
    s11_db = 20.0 * np.log10(np.abs(net_through.s[:, 0, 0]))
    s21_db = 20.0 * np.log10(np.abs(net_through.s[:, 1, 0]))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_through = OUT_DIR / "1D_SGBC_Through_skrf.csv"
    with out_through.open("w") as f:
        f.write("freq_hz,s11_db,s21_db\n")
        for freq, s11, s21 in zip(freqs, s11_db, s21_db):
            f.write(f"{freq:.12e},{s11:.12e},{s21:.12e}\n")
    print(f"Wrote {out_through} ({len(freqs)} frequencies, [{freqs[0]/1e9:.4f}, {freqs[-1]/1e9:.4f}] GHz)")

    # 2) PEC-backed lossy dielectric coating (eps_r=4, mu_r=1, sigma=10 S/m, d=5 mm)
    pec_layer = {
        "relative_permittivity": 4.0,
        "relative_permeability": 1.0,
        "bulk_conductivity": 10.0,
        "material_width": 0.005,
    }
    net_pec = build_multilayer_network(freqs, [pec_layer], pec_backed=True)
    s11_pec_db = 20.0 * np.log10(np.abs(net_pec.s[:, 0, 0]))

    out_pec = OUT_DIR / "1D_SGBC_Through_PEC_skrf.csv"
    with out_pec.open("w") as f:
        f.write("freq_hz,s11_db\n")
        for freq, s11 in zip(freqs, s11_pec_db):
            f.write(f"{freq:.12e},{s11:.12e}\n")
    print(f"Wrote {out_pec} ({len(freqs)} frequencies, [{freqs[0]/1e9:.4f}, {freqs[-1]/1e9:.4f}] GHz)")
    return out_through, out_pec


def generate_multiple_sgbc_references() -> tuple[Path, Path]:
    case_path = REPO / "testData" / "maxwellInputs" / "1D_SGBC_Multiple" / "1D_SGBC_Multiple.json"
    data = json.loads(case_path.read_text())
    freqs = probe_dft_frequencies(case_path, 0.2e9, 1.5e9)

    freq_obj = rf.Frequency.from_f(freqs, unit="Hz")
    air = Freespace(frequency=freq_obj)
    omega = 2.0 * np.pi * freq_obj.f

    def make_slab(eps_r: float, mu_r: float, sigma: float, width: float) -> rf.Network:
        eps_complex = eps_r - 1j * sigma / (omega * EPS0_SI)
        return Freespace(frequency=freq_obj, ep_r=eps_complex, mu_r=mu_r).line(d=width, unit="m")

    # 3 distinct SGBC boundaries at x = 1.00 m, x = 1.50 m, x = 2.00 m separated by 0.50 m vacuum gaps
    sgbc1 = make_slab(4.0, 1.0, 1.5, 0.003)
    sgbc2 = make_slab(2.5, 1.0, 0.5, 0.002) ** make_slab(5.0, 1.0, 2.5, 0.002)
    sgbc3 = make_slab(3.0, 1.0, 2.0, 0.0025)
    gap12 = air.line(d=0.50, unit="m")
    gap23 = air.line(d=0.50, unit="m")

    # 1) Passthrough cascade: air -> SGBC1 -> 0.5m air -> SGBC2 -> 0.5m air -> SGBC3 -> air
    net_pass = air.thru() ** sgbc1 ** gap12 ** sgbc2 ** gap23 ** sgbc3 ** air.thru()
    s11_db = 20.0 * np.log10(np.abs(net_pass.s[:, 0, 0]))
    s21_db = 20.0 * np.log10(np.abs(net_pass.s[:, 1, 0]))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_pass = OUT_DIR / "1D_SGBC_Multiple_skrf.csv"
    with out_pass.open("w") as f:
        f.write("freq_hz,s11_db,s21_db\n")
        for freq, s11, s21 in zip(freqs, s11_db, s21_db):
            f.write(f"{freq:.12e},{s11:.12e},{s21:.12e}\n")
    print(f"Wrote {out_pass} ({len(freqs)} frequencies, [{freqs[0]/1e9:.4f}, {freqs[-1]/1e9:.4f}] GHz)")

    # 2) PEC-backed cascade (SGBC3 right boundary set to PEC):
    net_pec = air.thru() ** sgbc1 ** gap12 ** sgbc2 ** gap23 ** sgbc3 ** air.short()
    s11_pec_db = 20.0 * np.log10(np.abs(net_pec.s[:, 0, 0]))

    out_pec = OUT_DIR / "1D_SGBC_Multiple_PEC_skrf.csv"
    with out_pec.open("w") as f:
        f.write("freq_hz,s11_db\n")
        for freq, s11 in zip(freqs, s11_pec_db):
            f.write(f"{freq:.12e},{s11:.12e}\n")
    print(f"Wrote {out_pec} ({len(freqs)} frequencies, [{freqs[0]/1e9:.4f}, {freqs[-1]/1e9:.4f}] GHz)")
    return out_pass, out_pec


if __name__ == "__main__":
    generate_pec_ram_reference()
    generate_passthrough_reference()
    generate_through_fast_references()
    generate_multiple_sgbc_references()
