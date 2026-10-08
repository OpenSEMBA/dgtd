#!/usr/bin/env python3
"""Standalone 1D nodal DG reference for a single-pole Lorentz slab.

This file does not call the C++ solver. It integrates

    dP/dt = J
    dJ/dt = -2 gamma J - omega_1^2 P + omega_p^2 E
    D = eps_inf E + P
    dD/dt = -dH/dx

with LGL nodes, an exact 1D Riemann flux, RK4, an SMA end at x=-5, a PEC
end at x=+5, and a total-field Gaussian at x=-3. Rates in a case JSON are
rad/s and are divided by c_SI, matching opensemba_dgtd. omega_1 = 0 is a
cold plasma.

Pass --sweep to write a few 1D cases, run the solver, and compare Ey.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C_SI = 299792458.0
REPO = Path(__file__).resolve().parents[1]
MESH_SRC = REPO / "testData/maxwellInputs/1D_Debye/1D_Debye.msh"
SOLVER = REPO / "build/gnu-release-mpi/bin/opensemba_dgtd"

# eps_inf, omega_p, omega_1, gamma in light-meters (1/length), folder, label
POLES = [
    (1.0, 2.0, 0.0, 0.4, "lorentz_plasma", r"plasma $\omega_p=2$, $\gamma=0.4$"),
    (2.0, 3.0, 2.0, 0.3, "lorentz_res2", r"$\omega_1=2$, $\omega_p=3$, $\gamma=0.3$"),
    (2.0, 2.5, 4.0, 0.1, "lorentz_res4", r"$\omega_1=4$, $\omega_p=2.5$, $\gamma=0.1$"),
    (1.0, 4.0, 0.0, 1.5, "lorentz_plasma_fast", r"plasma $\omega_p=4$, $\gamma=1.5$"),
]


def lgl_nodes_and_weights(order: int) -> tuple[np.ndarray, np.ndarray]:
    from numpy.polynomial.legendre import Legendre

    if order == 1:
        return np.array([-1.0, 1.0]), np.array([1.0, 1.0])
    basis = Legendre.basis(order)
    interior = basis.deriv().roots()
    nodes = np.concatenate(([-1.0], np.sort(interior.real), [1.0]))
    weights = 2.0 / (order * (order + 1) * basis(nodes) ** 2)
    return nodes, weights


def lgl_derivative_matrix(nodes: np.ndarray) -> np.ndarray:
    n = len(nodes)
    bary = np.ones(n)
    for j in range(n):
        diff = nodes[j] - nodes
        diff[j] = 1.0
        bary[j] = 1.0 / np.prod(diff)
    diff = nodes[:, None] - nodes[None, :]
    np.fill_diagonal(diff, 1.0)
    deriv = (bary[None, :] / bary[:, None]) / diff
    np.fill_diagonal(deriv, 0.0)
    np.fill_diagonal(deriv, -deriv.sum(axis=1))
    return deriv


def load_case(path: Path) -> dict:
    data = json.loads(path.read_text())
    opts = data["solver_options"]
    block = next(m for m in data["model"]["materials"] if "lorentz" in m)
    lorentz = block["lorentz"]
    source = data["sources"][0]
    mag = source["magnitude"]
    if mag.get("type", "gaussian") != "gaussian" or "frequency" in mag:
        raise SystemExit("This reference matches the unmodulated Gaussian planewave only.")
    inv_c = 1.0 / C_SI
    return {
        "path": path,
        "mesh": path.parent / data["model"]["filename"],
        "order": int(opts["order"]),
        "dt": float(opts["time_step"]),
        "t_end": float(opts["final_time"]),
        "save_every": float(data["probes"]["exporter"]["save_every"]),
        "eps_inf": float(lorentz["eps_inf"]),
        "omega_p": float(lorentz["omega_p"]) * inv_c,
        "omega_1": float(lorentz["omega_1"]) * inv_c,
        "gamma": float(lorentz["gamma"]) * inv_c,
        "spread": float(mag["spread"]),
        "tfsf_tag": int(source["tags"][0]),
        "lorentz_tags": [int(t) for t in block["tags"]],
    }


def load_mesh(path: Path, tfsf_tag: int, lorentz_tags: list[int]):
    nodes: dict[int, float] = {}
    lines: list[tuple[int, int, int]] = []
    tfsf_x = None
    in_nodes = False
    in_elems = False
    skip_count = False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line == "$Nodes":
            in_nodes = True
            skip_count = True
            continue
        if line == "$EndNodes":
            in_nodes = False
            continue
        if line == "$Elements":
            in_elems = True
            skip_count = True
            continue
        if line == "$EndElements":
            in_elems = False
            continue
        if skip_count:
            skip_count = False
            continue
        if in_nodes:
            parts = line.split()
            nodes[int(parts[0])] = float(parts[1])
            continue
        if not in_elems:
            continue
        parts = line.split()
        etype = int(parts[1])
        ntags = int(parts[2])
        tag = int(parts[3])
        conn = [int(v) for v in parts[3 + ntags :]]
        if etype == 15 and tag == tfsf_tag:
            tfsf_x = nodes[conn[0]]
        if etype == 1:
            lines.append((tag, conn[0], conn[1]))
    if tfsf_x is None:
        raise SystemExit(f"No TFSF point with tag {tfsf_tag} in {path}")
    elements = []
    for tag, n1, n2 in lines:
        xa, xb = nodes[n1], nodes[n2]
        if xa > xb:
            xa, xb = xb, xa
        elements.append((xa, xb, tag))
    elements.sort(key=lambda item: item[0])
    boundaries = np.array([elements[0][0]] + [xb for _, xb, _ in elements])
    material = np.array([1 if tag in lorentz_tags else 0 for _, _, tag in elements], dtype=int)
    tfsf_face = int(np.argmin(np.abs(boundaries[1:-1] - tfsf_x)))
    return boundaries, material, tfsf_face, float(tfsf_x)


def incident(t: float, x_face: float, mean: float, spread: float) -> tuple[float, float]:
    arg = (x_face - t) - mean
    ey = float(np.exp(-arg * arg / (2.0 * spread * spread)))
    return ey, ey


def riemann(e_l, h_l, eps_l, e_r, h_r, eps_r):
    z_l = 1.0 / np.sqrt(eps_l)
    z_r = 1.0 / np.sqrt(eps_r)
    a_l = e_l + z_l * h_l
    b_r = e_r - z_r * h_r
    denom = z_l + z_r
    e_star = (z_r * a_l + z_l * b_r) / denom
    h_star = (a_l - b_r) / denom
    return e_star, h_star


def rhs(d_field, h_field, p_field, j_field, t, geom):
    eps = geom["eps"]
    e_field = (d_field - p_field) / eps[:, None]
    inv_j = geom["inv_j"][:, None]
    deriv = geom["deriv"]
    d_dt = -(h_field @ deriv.T) * inv_j
    h_dt = -(e_field @ deriv.T) * inv_j
    p_dt = np.zeros_like(p_field)
    j_dt = np.zeros_like(j_field)
    mask = geom["lorentz"]
    if np.any(mask):
        p_dt[mask] = j_field[mask]
        j_dt[mask] = (
            -2.0 * geom["gamma"][mask, None] * j_field[mask]
            - geom["omega_1_sq"][mask, None] * p_field[mask]
            + geom["omega_p_sq"][mask, None] * e_field[mask]
        )

    e_l, e_r = e_field[:-1, -1], e_field[1:, 0]
    h_l, h_r = h_field[:-1, -1], h_field[1:, 0]
    e_star, h_star = riemann(e_l, h_l, eps[:-1], e_r, h_r, eps[1:])
    e_star_r = e_star.copy()
    h_star_r = h_star.copy()
    face = geom["tfsf_face"]
    e_inc, h_inc = incident(t, geom["tfsf_x"], geom["mean"], geom["spread"])
    e_sc, h_sc = riemann(
        e_l[face], h_l[face], 1.0,
        e_r[face] - e_inc, h_r[face] - h_inc, 1.0,
    )
    e_star_l = e_star.copy()
    h_star_l = h_star.copy()
    e_star_l[face] = e_sc
    h_star_l[face] = h_sc
    e_star_r[face] = e_sc + e_inc
    h_star_r[face] = h_sc + h_inc

    w = geom["weights"]
    jac = geom["j"]
    d_dt[1:, 0] += (h_star_r - h_r) / jac[1:] / w[0]
    h_dt[1:, 0] += (e_star_r - e_r) / jac[1:] / w[0]
    d_dt[:-1, -1] -= (h_star_l - h_l) / jac[:-1] / w[-1]
    h_dt[:-1, -1] -= (e_star_l - e_l) / jac[:-1] / w[-1]

    e_int, h_int = e_field[0, 0], h_field[0, 0]
    b_out = e_int - h_int
    d_dt[0, 0] += (-0.5 * b_out - h_int) / jac[0] / w[0]
    h_dt[0, 0] += (0.5 * b_out - e_int) / jac[0] / w[0]

    e_int = e_field[-1, -1]
    h_dt[-1, -1] -= (0.0 - e_int) / jac[-1] / w[-1]
    return d_dt, h_dt, p_dt, j_dt


def rk4(state, t, dt, geom):
    stages = []
    times = (t, t + 0.5 * dt, t + 0.5 * dt, t + dt)
    scales = (0.0, 0.5, 0.5, 1.0)
    for scale, time in zip(scales, times):
        if scale == 0.0:
            packed = state
        else:
            prev = stages[-1]
            packed = tuple(s + scale * dt * k for s, k in zip(state, prev))
        stages.append(rhs(*packed, time, geom))
    k1, k2, k3, k4 = stages

    def mix(s, *ks):
        return s + (dt / 6.0) * (ks[0] + 2.0 * ks[1] + 2.0 * ks[2] + ks[3])

    return tuple(
        mix(s, k1[i], k2[i], k3[i], k4[i]) for i, s in enumerate(state)
    )


def build_geom(case: dict):
    boundaries, material, tfsf_face, tfsf_x = load_mesh(
        case["mesh"], case["tfsf_tag"], case["lorentz_tags"]
    )
    nodes, weights = lgl_nodes_and_weights(case["order"])
    deriv = lgl_derivative_matrix(nodes)
    jac = 0.5 * np.diff(boundaries)
    xa = boundaries[:-1]
    x_nodes = 0.5 * (xa + boundaries[1:])[:, None] + jac[:, None] * nodes[None, :]
    eps = np.where(material == 1, case["eps_inf"], 1.0)
    lorentz = material == 1
    return {
        "boundaries": boundaries,
        "x_nodes": x_nodes,
        "j": jac,
        "inv_j": 1.0 / jac,
        "weights": weights,
        "deriv": deriv,
        "eps": eps,
        "lorentz": lorentz,
        "omega_p_sq": np.where(lorentz, case["omega_p"] ** 2, 0.0),
        "omega_1_sq": np.where(lorentz, case["omega_1"] ** 2, 0.0),
        "gamma": np.where(lorentz, case["gamma"], 0.0),
        "tfsf_face": tfsf_face,
        "tfsf_x": tfsf_x,
        "mean": tfsf_x - 5.0 * case["spread"] * np.sqrt(2.0),
        "spread": case["spread"],
    }


def run_reference(case: dict, geom: dict):
    n_steps = int(round(case["t_end"] / case["dt"]))
    dt = case["t_end"] / n_steps
    save_every_steps = max(1, int(round(case["save_every"] / dt)))
    ne, n_p = geom["x_nodes"].shape
    state = tuple(np.zeros((ne, n_p)) for _ in range(4))
    frames = []
    times = []
    t = 0.0
    for step in range(n_steps + 1):
        if step % save_every_steps == 0 or step == n_steps:
            d_field, _, p_field, _ = state
            frames.append(((d_field - p_field) / geom["eps"][:, None]).copy())
            times.append(t)
        if step == n_steps:
            break
        state = rk4(state, t, dt, geom)
        t += dt
    return np.asarray(times), frames, dt, n_steps


def analytic_vacuum(x, t, mean, spread):
    arg = (x - t) - mean
    return np.exp(-arg * arg / (2.0 * spread * spread))


def vacuum_check(geom, times, frames, spread):
    idx = int(np.argmin(np.abs(times - 8.0)))
    x = geom["x_nodes"]
    window = (x > geom["tfsf_x"] + 0.3) & (x < 2.5)
    exact = analytic_vacuum(x[window], times[idx], geom["mean"], spread)
    num = frames[idx][window]
    rel = np.linalg.norm(num - exact) / max(np.linalg.norm(exact), 1e-30)
    return float(times[idx]), float(rel)


def pvd_datasets(pvd: Path) -> list[tuple[float, Path]]:
    root = ET.parse(pvd).getroot()
    return [
        (float(ds.attrib["timestep"]), pvd.parent / ds.attrib["file"])
        for ds in root.iter("DataSet")
    ]


def read_ey(vtu: Path) -> tuple[np.ndarray, np.ndarray]:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLPUnstructuredGridReader()
    reader.SetFileName(str(vtu))
    reader.Update()
    grid = reader.GetOutput()
    pts = vtk_to_numpy(grid.GetPoints().GetData())
    field = grid.GetPointData().GetArray("E")
    if field is None:
        raise SystemExit(f"{vtu} has no E array")
    return pts[:, 0], vtk_to_numpy(field)[:, 1]


def sample_reference(x_query, x_nodes, ey, boundaries):
    values = np.empty_like(x_query)
    for i, x in enumerate(x_query):
        e = int(np.searchsorted(boundaries[1:], x, side="right"))
        e = min(max(e, 0), len(boundaries) - 2)
        order = np.argsort(x_nodes[e])
        xs = x_nodes[e, order]
        ys = ey[e, order]
        if x <= xs[0]:
            values[i] = ys[0]
        elif x >= xs[-1]:
            values[i] = ys[-1]
        else:
            values[i] = np.interp(x, xs, ys)
    return values


def compare(case, geom, times, frames, pvd: Path, out_png: Path):
    sets = pvd_datasets(pvd)
    rows = []
    set_times = np.array([t for t, _ in sets])
    for t_ref, frame in zip(times, frames):
        j = int(np.argmin(np.abs(set_times - t_ref)))
        t_dg, vtu = sets[j]
        if abs(t_dg - t_ref) > 0.51 * case["save_every"]:
            continue
        x_dg, ey_dg = read_ey(vtu)
        ey_ref = sample_reference(x_dg, geom["x_nodes"], frame, geom["boundaries"])
        den = np.linalg.norm(ey_ref)
        rel = np.linalg.norm(ey_dg - ey_ref) / den if den > 1e-8 else np.nan
        rows.append((t_ref, rel, x_dg, ey_dg, ey_ref))
    usable = [r for r in rows if np.isfinite(r[1])]
    if not usable:
        print("No overlapping export times with a nonzero reference field.")
        return []
    peak = max(np.linalg.norm(r[3]) for r in usable)
    developed = [r for r in usable if r[0] >= 6.0 and np.linalg.norm(r[3]) > 0.05 * peak]
    worst = max(developed or usable, key=lambda r: r[1])
    print(f"Compared {len(usable)} frames against {pvd}")
    print(f"Largest relative L2 for t>=6 is {worst[1]:.3e} at t={worst[0]:.3f}")
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(worst[2], worst[4], label="reference")
    ax.plot(worst[2], worst[3], linestyle="--", label="dgtd")
    ax.axvline(geom["tfsf_x"], color="0.4", linestyle=":")
    ax.axvline(3.0, color="0.2", linestyle="--")
    ax.set_xlabel("x")
    ax.set_ylabel("Ey")
    ax.set_title(f"t = {worst[0]:.2f}, relative L2 = {worst[1]:.3e}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    print(f"Wrote {out_png}")
    return usable


def write_case(name, eps_inf, omega_p, omega_1, gamma):
    folder = REPO / "testData/maxwellInputs" / name
    folder.mkdir(parents=True, exist_ok=True)
    shutil.copy(MESH_SRC, folder / f"{name}.msh")
    data = {
        "solver_options": {
            "upwind_alpha": 1.0,
            "final_time": 25.0,
            "order": 4,
            "time_step": 0.01,
        },
        "model": {
            "filename": f"{name}.msh",
            "materials": [
                {"tags": [1, 2], "type": "vacuum"},
                {
                    "tags": [3],
                    "lorentz": {
                        "eps_inf": eps_inf,
                        "omega_p": omega_p * C_SI,
                        "omega_1": omega_1 * C_SI,
                        "gamma": gamma * C_SI,
                    },
                },
            ],
            "boundaries": [
                {"tags": [1], "type": "SMA"},
                {"tags": [2], "type": "PEC"},
            ],
        },
        "probes": {"exporter": {"save_every": 1.0}},
        "sources": [
            {
                "type": "planewave",
                "polarization": [0.0, 1.0, 0.0],
                "propagation": [1.0, 0.0, 0.0],
                "tags": [3],
                "magnitude": {"type": "gaussian", "spread": 0.6},
            }
        ],
    }
    path = folder / f"{name}.json"
    path.write_text(json.dumps(data, indent=2) + "\n")
    return path


def band(rows, t0, t1):
    group = [r[1] for r in rows if t0 <= r[0] <= t1 and np.isfinite(r[1])]
    return None if not group else float(max(group))


def sweep():
    if not SOLVER.is_file():
        raise SystemExit(f"Solver not found: {SOLVER}")
    if not MESH_SRC.is_file():
        raise SystemExit(f"Mesh not found: {MESH_SRC}")
    summary = []
    panels = []
    for eps_inf, omega_p, omega_1, gamma, name, label in POLES:
        case_path = write_case(name, eps_inf, omega_p, omega_1, gamma)
        print(f"\n=== {name} ===")
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
        early = band(rows, 5.0, 9.0)
        late = band(rows, 12.0, 20.0)
        snap = min(rows, key=lambda r: abs(r[0] - 14.0)) if rows else None
        print(
            f"  steps {n_steps}, analytic t={t_check:.2f} rel {analytic:.3e}, "
            f"vacuum max {early}, slab max {late}"
        )
        summary.append((label, early, late, analytic))
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
    fig.suptitle("Lorentz slab, same Gaussian, rates in light-meters")
    fig.tight_layout()
    out = REPO / "exports/SimulationData/single-core/lorentz_1d/lorentz_pole_sweep.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"\nWrote {out}")
    for label, early, late, analytic in summary:
        ev = "n/a" if early is None else f"{early:.3e}"
        lv = "n/a" if late is None else f"{late:.3e}"
        print(f"{label}: vacuum {ev}, slab {lv}, analytic {analytic:.3e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", nargs="?", type=Path)
    parser.add_argument("--sweep", action="store_true")
    parser.add_argument("--pvd", type=Path)
    parser.add_argument("--png", type=Path)
    args = parser.parse_args()
    if args.sweep:
        sweep()
        return
    if args.case is None:
        raise SystemExit("Pass a case JSON, or --sweep.")
    case = load_case(args.case)
    geom = build_geom(case)
    print("1D Lorentz reference")
    print(
        f"  eps_inf {case['eps_inf']}, omega_p {case['omega_p']:.6g}, "
        f"omega_1 {case['omega_1']:.6g}, gamma {case['gamma']:.6g} (light-meters)"
    )
    times, frames, dt, n_steps = run_reference(case, geom)
    t_check, rel = vacuum_check(geom, times, frames, case["spread"])
    print(f"  {n_steps} steps, dt {dt}, analytic t={t_check:.3f} rel {rel:.3e}")
    if args.pvd and args.pvd.is_file():
        png = args.png or (args.pvd.parent / "compare.png")
        compare(case, geom, times, frames, args.pvd, png)


if __name__ == "__main__":
    main()
