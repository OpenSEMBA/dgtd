#!/usr/bin/env python3
"""Yuan 2010 conductor-backed plasma slab, run with opensemba_dgtd.

Reproduces the 1D case in lmdiazangulo/PyDG1D (branch lorentz-model),
test/test_lorentz_plasma_slab.py and examples/plot_yuan2010_dgtd_10cm.py.

C. X. Yuan, Z. X. Zhou and H. G. Sun, "Reflection Properties of
Electromagnetic Wave in a Bounded Plasma Slab", IEEE Trans. Plasma Sci.,
vol. 38, no. 12, pp. 3348-3355, 2010. Cold plasma, f_p = 2 GHz,
nu_en = 5e9 s^-1, slab thickness 10 cm, PEC backing. The reflected power
is their Eqs. (10)-(12). The time-domain model is the cold-plasma limit of
the single-pole Lorentz ODE (omega_1 = 0, eps_inf = 1, gamma = nu_en / 2).

Geometry matches that script: order 2, dx = 5 mm, a right-going Gaussian of
width 0.02 m centered at x = 0.25 m, probe at x = 0.6 m, plasma on
[0.9, 1.0] m, PEC at x = 1 m. A second vacuum run on [0, 3] m supplies the
incident trace. The left wall is SMA, which in 1D is the same outgoing
characteristic condition as their ABC.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C_SI = 299792458.0
REPO = Path(__file__).resolve().parents[1]
SOLVER = REPO / "build/gnu-release-mpi/bin/opensemba_dgtd"

PLASMA_FREQUENCY = 2.0e9
COLLISION_FREQUENCY = 5.0e9
SLAB_THICKNESS = 0.10
FREQUENCIES = np.linspace(2.0e9, 4.0e9, 201)
XMAX_PLASMA = 1.0
XMAX_VACUUM = 3.0
N_PLASMA = 200
DX = XMAX_PLASMA / N_PLASMA
PROBE_X = 0.6
SOURCE_X = 0.25
SOURCE_WIDTH = 0.02
T_END = 2.5
DT = 5.0e-4
ORDER = 2


def write_mesh(path: Path, xmax: float, slab_from: float | None, right_tag: int):
    n_elem = int(round(xmax / DX))
    xs = [i * DX for i in range(n_elem + 1)]
    nodes = [(i + 1, xs[i]) for i in range(n_elem + 1)]
    lines = []
    eid = 1
    lines.append(f"{eid} 15 2 1 1 1")
    eid += 1
    lines.append(f"{eid} 15 2 {right_tag} {right_tag} {n_elem + 1}")
    eid += 1
    for e in range(n_elem):
        x_left = xs[e]
        tag = 2 if slab_from is not None and x_left >= slab_from - 1e-12 else 1
        lines.append(f"{eid} 1 2 {tag} {tag} {e + 1} {e + 2}")
        eid += 1
    text = [
        "$MeshFormat",
        "2.2 0 8",
        "$EndMeshFormat",
        "$PhysicalNames",
        "1",
        '1 1 "vacuum"',
        "$EndPhysicalNames",
        "$Nodes",
        str(len(nodes)),
    ]
    text.extend(f"{i} {x:.16g} 0 0" for i, x in nodes)
    text.append("$EndNodes")
    text.append("$Elements")
    text.append(str(len(lines)))
    text.extend(lines)
    text.append("$EndElements")
    path.write_text("\n".join(text) + "\n")


def write_case(folder: Path, name: str, xmax: float, plasma: bool, right: str):
    folder.mkdir(parents=True, exist_ok=True)
    write_mesh(
        folder / f"{name}.msh",
        xmax,
        XMAX_PLASMA - SLAB_THICKNESS if plasma else None,
        2,
    )
    materials = [{"tags": [1], "type": "vacuum"}]
    if plasma:
        materials.append(
            {
                "tags": [2],
                "lorentz": {
                    "eps_inf": 1.0,
                    "omega_p": 2.0 * np.pi * PLASMA_FREQUENCY,
                    "omega_1": 0.0,
                    "gamma": 0.5 * COLLISION_FREQUENCY,
                },
            }
        )
    data = {
        "solver_options": {
            "evolution_operator": "global",
            "upwind_alpha": 1.0,
            "final_time": T_END,
            "order": ORDER,
            "time_step": DT,
        },
        "model": {
            "filename": f"{name}.msh",
            "materials": materials,
            "boundaries": [
                {"tags": [1], "type": "SMA"},
                {"tags": [2], "type": right},
            ],
        },
        "probes": {"point": [{"position": [PROBE_X], "steps": 1}]},
        "sources": [
            {
                "type": "initial",
                "field_type": "electric",
                "center": [SOURCE_X],
                "polarization": [0.0, 1.0, 0.0],
                "dimension": 1,
                "magnitude": {"type": "gaussian", "spread": SOURCE_WIDTH},
            },
            {
                "type": "initial",
                "field_type": "magnetic",
                "center": [SOURCE_X],
                "polarization": [0.0, 0.0, 1.0],
                "dimension": 1,
                "magnitude": {"type": "gaussian", "spread": SOURCE_WIDTH},
            },
        ],
    }
    path = folder / f"{name}.json"
    path.write_text(json.dumps(data, indent=2) + "\n")
    return path


def analytical_reflection_db(frequencies: np.ndarray) -> np.ndarray:
    """Yuan 2010 Eqs. (10)-(12), reflected power in dB."""
    omega = 2.0 * np.pi * frequencies
    omega_p = 2.0 * np.pi * PLASMA_FREQUENCY
    eps = 1.0 - omega_p**2 / (omega * (omega - 1j * COLLISION_FREQUENCY))
    root = np.sqrt(eps)
    # Passive branch: decaying wave into the slab, Re(sqrt(eps)) >= 0.
    root = np.where(np.real(root) < 0.0, -root, root)
    arg = 1j * omega * SLAB_THICKNESS * root / C_SI
    th = np.tanh(arg)
    reflection = (th - root) / (th + root)
    return 10.0 * np.log10(np.abs(reflection) ** 2)


def run_solver(case: Path):
    log = REPO / "exports/SimulationData/single-core" / case.stem / "solver.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as out:
        subprocess.run(
            ["mpirun", "-np", "1", str(SOLVER), "-i", str(case)],
            cwd=REPO,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return log


def load_ey(case_name: str):
    path = (
        REPO
        / "exports/SimulationData/single-core"
        / case_name
        / "PointProbes/PointProbe0.dat"
    )
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 7:
            continue
        try:
            values = [float(v) for v in parts]
        except ValueError:
            continue
        rows.append(values)
    data = np.asarray(rows)
    if data.size == 0:
        raise SystemExit(f"No samples in {path}")
    return data[:, 0], data[:, 2]


def reflection_db(time, total, incident, frequencies):
    dt = float(np.median(np.diff(time)))

    def dft(values):
        phase = np.exp(-2j * np.pi * frequencies[:, None] * time[None, :])
        return (phase * values[None, :]).sum(axis=1) * dt

    reflected = dft(total - incident)
    inc = dft(incident)
    return 20.0 * np.log10(np.abs(reflected) / np.abs(inc))


def main():
    if not SOLVER.is_file():
        raise SystemExit(f"Solver not found: {SOLVER}")
    vacuum = write_case(
        REPO / "testData/maxwellInputs/yuan2010_vacuum",
        "yuan2010_vacuum",
        XMAX_VACUUM,
        plasma=False,
        right="SMA",
    )
    plasma = write_case(
        REPO / "testData/maxwellInputs/yuan2010_plasma",
        "yuan2010_plasma",
        XMAX_PLASMA,
        plasma=True,
        right="PEC",
    )
    print("Running vacuum reference")
    run_solver(vacuum)
    print("Running plasma slab")
    run_solver(plasma)

    t_vac, ey_vac = load_ey("yuan2010_vacuum")
    t_pla, ey_pla = load_ey("yuan2010_plasma")
    n = min(len(t_vac), len(t_pla))
    if np.max(np.abs(t_vac[:n] - t_pla[:n])) > 1e-18:
        raise SystemExit("Probe time bases differ.")
    time, incident, total = t_vac[:n], ey_vac[:n], ey_pla[:n]
    peak = time[int(np.argmax(np.abs(incident)))]
    print(
        f"Incident Ey peak {incident[int(np.argmax(np.abs(incident)))]:.4f} "
        f"at t = {peak * 1e9:.3f} ns "
        f"(expected near {(PROBE_X - SOURCE_X) / C_SI * 1e9:.3f} ns)"
    )

    expected = analytical_reflection_db(FREQUENCIES)
    numerical = reflection_db(time, total, incident, FREQUENCIES)
    error = numerical - expected
    print(
        f"dgtd vs Yuan2010: max |error| = {np.max(np.abs(error)):.3f} dB, "
        f"RMS = {np.sqrt(np.mean(error**2)):.3f} dB"
    )
    print(
        f"analytical dip {expected.min():.2f} dB "
        f"at {FREQUENCIES[np.argmin(expected)] / 1e9:.3f} GHz"
    )
    print(
        f"dgtd dip {numerical.min():.2f} dB "
        f"at {FREQUENCIES[np.argmin(numerical)] / 1e9:.3f} GHz"
    )

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(FREQUENCIES / 1e9, expected, color="black", linewidth=2.0,
            label="Yuan 2010, Eqs. (10)-(12)")
    ax.plot(FREQUENCIES / 1e9, numerical, color="tab:red", linewidth=1.4,
            linestyle="--", label="opensemba_dgtd")
    ax.set_xlim(2.0, 4.0)
    ax.set_ylim(-30.0, 0.0)
    ax.set_xlabel("Frequency / GHz")
    ax.set_ylabel("Total reflection [dB]")
    ax.set_title("10 cm plasma slab over a conductor")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", framealpha=0.95)
    fig.tight_layout()
    out = REPO / "exports/SimulationData/single-core/yuan2010_plasma/yuan2010_reflection.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
