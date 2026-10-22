#!/usr/bin/env python3
"""Run several Debye poles on the 1D slab and compare each with the reference.

JSON tau is seconds. The value written here is tau_solver / c_SI, so the
numbers below are the poles the time stepper actually integrates, in
light-meters.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from debye_1d_reference import (
    C_SI,
    REPO,
    build_geom,
    compare,
    load_case,
    run_reference,
    vacuum_check,
)

SOLVER = REPO / "build/gnu-release-mpi/bin/opensemba_dgtd"
MESH_SRC = REPO / "testData/maxwellInputs/1D_Debye/1D_Debye.msh"
TEMPLATE = json.loads(
    (REPO / "testData/maxwellInputs/1D_Debye/1D_Debye.json").read_text()
)

# eps_inf, eps_s, tau in light-meters, short label
POLES = [
    (2.0, 30.0, 0.2, "eps30_tau0p2", r"$\varepsilon_s=30$, $\tau=0.2$"),
    (2.0, 30.0, 0.8, "eps30_tau0p8", r"$\varepsilon_s=30$, $\tau=0.8$"),
    (2.0, 30.0, 4.0, "eps30_tau4", r"$\varepsilon_s=30$, $\tau=4$"),
    (2.0, 6.0, 0.8, "eps6_tau0p8", r"$\varepsilon_s=6$, $\tau=0.8$"),
]

T_END = 25.0
DT = 0.01
SAVE_EVERY = 1.0
ORDER = 4
SNAPSHOT = 14.0


def write_case(name: str, eps_inf: float, eps_s: float, tau_solver: float) -> Path:
    folder = REPO / "testData/maxwellInputs" / name
    folder.mkdir(parents=True, exist_ok=True)
    shutil.copy(MESH_SRC, folder / f"{name}.msh")
    data = json.loads(json.dumps(TEMPLATE))
    data["solver_options"]["final_time"] = T_END
    data["solver_options"]["time_step"] = DT
    data["solver_options"]["order"] = ORDER
    data["model"]["filename"] = f"{name}.msh"
    data["probes"]["exporter"]["save_every"] = SAVE_EVERY
    data["model"]["materials"][1]["debye"] = {
        "eps_inf": eps_inf,
        "eps_s": eps_s,
        "tau": tau_solver / C_SI,
    }
    path = folder / f"{name}.json"
    path.write_text(json.dumps(data, indent=2) + "\n")
    return path


def band(rows, t0, t1) -> float | None:
    group = [r[1] for r in rows if t0 <= r[0] <= t1 and np.isfinite(r[1])]
    if not group:
        return None
    return float(max(group))


def main():
    if not SOLVER.is_file():
        raise SystemExit(f"Solver not found: {SOLVER}")
    summary = []
    panels = []
    for eps_inf, eps_s, tau_solver, name, label in POLES:
        case_path = write_case(name, eps_inf, eps_s, tau_solver)
        print(f"\n=== {name}  eps_s={eps_s}  tau={tau_solver} light-meters ===")
        log = REPO / "exports/SimulationData/single-core" / name / "solver.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as out:
            subprocess.run(
                ["mpirun", "-np", "1", str(SOLVER), "-i", str(case_path)],
                cwd=REPO,
                stdout=out,
                stderr=subprocess.STDOUT,
                check=True,
            )
        case = load_case(case_path)
        geom = build_geom(case)
        times, frames, _, n_steps = run_reference(case, geom)
        t_check, analytic = vacuum_check(geom, times, frames, case["spread"])
        pvd = REPO / "exports/ParaView/single-core" / f"{name}.msh" / f"{name}.msh.pvd"
        png = REPO / "exports/SimulationData/single-core" / name / "compare.png"
        rows = compare(case, geom, times, frames, pvd, png) or []
        late = band(rows, 12.0, 20.0)
        early = band(rows, 5.0, 9.0)
        snap = min(rows, key=lambda r: abs(r[0] - SNAPSHOT)) if rows else None
        print(
            f"  steps {n_steps}, analytic t={t_check:.2f} rel {analytic:.3e}, "
            f"vacuum max rel {early}, slab max rel {late}"
        )
        summary.append((label, tau_solver, eps_s, early, late, analytic))
        if snap is not None:
            panels.append((label, snap))

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
    for ax, (label, snap) in zip(axes.ravel(), panels):
        t_ref, rel, x_dg, ey_dg, ey_ref = snap
        ax.plot(x_dg, ey_ref, label="reference")
        ax.plot(x_dg, ey_dg, linestyle="--", label="dgtd")
        ax.axvline(-3.0, color="0.5", linestyle=":", linewidth=1)
        ax.axvline(3.0, color="0.2", linestyle="--", linewidth=1)
        ax.set_title(f"{label}\n$t={t_ref:.0f}$, rel $L_2={rel:.2e}$")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    for ax in axes[1, :]:
        ax.set_xlabel("x")
    for ax in axes[:, 0]:
        ax.set_ylabel("Ey")
    fig.suptitle(r"Debye slab, $\varepsilon_\infty=2$, same pulse")
    fig.tight_layout()
    out = REPO / "exports/SimulationData/single-core/1D_Debye/debye_pole_sweep.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"\nWrote {out}")
    print(f"{'pole':<28} {'vacuum max':>12} {'slab max':>12}")
    for label, tau, eps_s, early, late, analytic in summary:
        ev = "n/a" if early is None else f"{early:.3e}"
        lv = "n/a" if late is None else f"{late:.3e}"
        print(f"eps_s={eps_s:<4} tau={tau:<4} {ev:>12} {lv:>12}")


if __name__ == "__main__":
    main()
