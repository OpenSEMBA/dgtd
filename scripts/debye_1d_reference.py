#!/usr/bin/env python3
"""Fast 1D nodal DG reference for testData/maxwellInputs/1D_Debye.

The loop in the original script walked every element in Python on every
Runge-Kutta stage. This version keeps that scheme (LGL nodes, strong-form
DG, exact 1D Riemann flux, RK4) and applies it with array operations.

Parameters are read from the case JSON and interpreted the way opensemba_dgtd
does:

- solver time is light-meters (c = 1)
- debye.tau is seconds, so the pole used here is tau * c_SI
- the planewave is an unmodulated Gaussian; spread is the solver sigma
- the pulse is placed 5*sigma*sqrt(2) upstream of the total-field face

The face law here is the exact impedance Riemann flux. The solver's upwind
flux does not use Z(eps) on a material jump, so a Debye-interface reflection
is not expected to match digit for digit.
"""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C_SI = 299792458.0
REPO = Path(__file__).resolve().parents[1]
DEFAULT_CASE = REPO / "testData" / "maxwellInputs" / "1D_Debye" / "1D_Debye.json"


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
    debye = next(m["debye"] for m in data["model"]["materials"] if "debye" in m)
    source = data["sources"][0]
    mag = source["magnitude"]
    if mag.get("type", "gaussian") != "gaussian" or "frequency" in mag:
        raise SystemExit("This reference matches the unmodulated Gaussian planewave only.")
    tau_json = float(debye["tau"])
    return {
        "path": path,
        "mesh": path.parent / data["model"]["filename"],
        "order": int(opts["order"]),
        "dt": float(opts["time_step"]),
        "t_end": float(opts["final_time"]),
        "save_every": float(data["probes"]["exporter"]["save_every"]),
        "eps_inf": float(debye["eps_inf"]),
        "eps_s": float(debye["eps_s"]),
        "tau_json": tau_json,
        "tau": tau_json * C_SI,
        "spread": float(mag["spread"]),
        "tfsf_tag": int(source["tags"][0]),
        "debye_tags": [int(t) for t in next(
            m["tags"] for m in data["model"]["materials"] if "debye" in m
        )],
    }


def load_mesh(path: Path, tfsf_tag: int, debye_tags: list[int]):
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
    material = np.array([1 if tag in debye_tags else 0 for _, _, tag in elements], dtype=int)
    faces = 0.5 * (boundaries[:-1] + boundaries[1:])
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


def rhs(d_field, h_field, p_field, t, geom):
    eps = geom["eps"]
    e_field = (d_field - p_field) / eps[:, None]
    inv_j = geom["inv_j"][:, None]
    deriv = geom["deriv"]
    d_dt = -(h_field @ deriv.T) * inv_j
    h_dt = -(e_field @ deriv.T) * inv_j
    p_dt = np.zeros_like(p_field)
    debye = geom["debye"]
    if np.any(debye):
        p_dt[debye] = geom["eps_d"] / geom["tau"] * e_field[debye] - p_field[debye] / geom["tau"]

    e_l, e_r = e_field[:-1, -1], e_field[1:, 0]
    h_l, h_r = h_field[:-1, -1], h_field[1:, 0]
    eps_l, eps_r = eps[:-1], eps[1:]
    e_star, h_star = riemann(e_l, h_l, eps_l, e_r, h_r, eps_r)
    e_star_r, h_star_r = e_star.copy(), h_star.copy()

    # Scattered Riemann on the total-field face. The right element stores the
    # total field, so its star state carries the incident wave as well.
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
    j = geom["j"]
    d_dt[1:, 0] += (h_star_r - h_r) / j[1:] / w[0]
    h_dt[1:, 0] += (e_star_r - e_r) / j[1:] / w[0]
    d_dt[:-1, -1] -= (h_star_l - h_l) / j[:-1] / w[-1]
    h_dt[:-1, -1] -= (e_star_l - e_l) / j[:-1] / w[-1]

    e_int, h_int = e_field[0, 0], h_field[0, 0]
    b_out = e_int - h_int
    e_sma, h_sma = 0.5 * b_out, -0.5 * b_out
    d_dt[0, 0] += (h_sma - h_int) / j[0] / w[0]
    h_dt[0, 0] += (e_sma - e_int) / j[0] / w[0]

    e_int, h_int = e_field[-1, -1], h_field[-1, -1]
    d_dt[-1, -1] -= (h_int - h_int) / j[-1] / w[-1]
    h_dt[-1, -1] -= (0.0 - e_int) / j[-1] / w[-1]
    return d_dt, h_dt, p_dt


def rk4(d_field, h_field, p_field, t, dt, geom):
    stages = []
    state = (d_field, h_field, p_field)
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
    return (
        mix(d_field, k1[0], k2[0], k3[0], k4[0]),
        mix(h_field, k1[1], k2[1], k3[1], k4[1]),
        mix(p_field, k1[2], k2[2], k3[2], k4[2]),
    )


def build_geom(case: dict):
    boundaries, material, tfsf_face, tfsf_x = load_mesh(
        case["mesh"], case["tfsf_tag"], case["debye_tags"]
    )
    nodes, weights = lgl_nodes_and_weights(case["order"])
    deriv = lgl_derivative_matrix(nodes)
    j = 0.5 * np.diff(boundaries)
    xa = boundaries[:-1]
    x_nodes = 0.5 * (xa + boundaries[1:])[:, None] + j[:, None] * nodes[None, :]
    eps = np.where(material == 1, case["eps_inf"], 1.0)
    mean = tfsf_x - 5.0 * case["spread"] * np.sqrt(2.0)
    return {
        "boundaries": boundaries,
        "x_nodes": x_nodes,
        "j": j,
        "inv_j": 1.0 / j,
        "weights": weights,
        "deriv": deriv,
        "eps": eps,
        "debye": material == 1,
        "eps_d": case["eps_s"] - case["eps_inf"],
        "tau": case["tau"],
        "tfsf_face": tfsf_face,
        "tfsf_x": tfsf_x,
        "mean": mean,
        "spread": case["spread"],
        "material": material,
    }


def run_reference(case: dict, geom: dict):
    n_steps = int(round(case["t_end"] / case["dt"]))
    dt = case["t_end"] / n_steps
    save_every_steps = max(1, int(round(case["save_every"] / dt)))
    ne, n_p = geom["x_nodes"].shape
    d_field = np.zeros((ne, n_p))
    h_field = np.zeros((ne, n_p))
    p_field = np.zeros((ne, n_p))
    frames = []
    times = []
    t = 0.0
    for step in range(n_steps + 1):
        if step % save_every_steps == 0 or step == n_steps:
            frames.append(((d_field - p_field) / geom["eps"][:, None]).copy())
            times.append(t)
        if step == n_steps:
            break
        d_field, h_field, p_field = rk4(d_field, h_field, p_field, t, dt, geom)
        t += dt
    return np.asarray(times), frames, dt, n_steps


def analytic_vacuum(x, t, mean, spread):
    arg = (x - t) - mean
    return np.exp(-arg * arg / (2.0 * spread * spread))


def vacuum_check(geom, times, frames, spread):
    target = 8.0
    idx = int(np.argmin(np.abs(times - target)))
    x = geom["x_nodes"]
    ey = frames[idx]
    window = (x > geom["tfsf_x"] + 0.3) & (x < 2.5)
    exact = analytic_vacuum(x[window], times[idx], geom["mean"], spread)
    num = ey[window]
    rel = np.linalg.norm(num - exact) / max(np.linalg.norm(exact), 1e-30)
    return times[idx], rel


def pvd_datasets(pvd: Path) -> list[tuple[float, Path]]:
    root = ET.parse(pvd).getroot()
    out = []
    for ds in root.iter("DataSet"):
        out.append((float(ds.attrib["timestep"]), pvd.parent / ds.attrib["file"]))
    return out


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
    ey = vtk_to_numpy(field)[:, 1]
    return pts[:, 0], ey


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
    for t_ref, frame in zip(times, frames):
        j = int(np.argmin(abs(np.array([t for t, _ in sets]) - t_ref)))
        t_dg, vtu = sets[j]
        if abs(t_dg - t_ref) > 0.51 * case["save_every"]:
            continue
        x_dg, ey_dg = read_ey(vtu)
        ey_ref = sample_reference(x_dg, geom["x_nodes"], frame, geom["boundaries"])
        num = np.linalg.norm(ey_dg)
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
    early = [r for r in usable if 5.0 <= r[0] <= 9.0]
    late = [r for r in usable if 12.0 <= r[0] <= 20.0]
    def band(name, group):
        if not group:
            print(f"{name}: no snapshot in range")
            return
        errs = [r[1] for r in group]
        print(f"{name}: relative L2 {min(errs):.3e} to {max(errs):.3e} over {len(group)} frames")
    print(f"Compared {len(usable)} frames against {pvd}")
    band("Vacuum transit t=5..9", early)
    band("After the slab t=12..20", late)
    print(f"Largest relative L2 is {worst[1]:.3e} at t={worst[0]:.3f}")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(worst[2], worst[4], label="reference")
    ax.plot(worst[2], worst[3], label="dgtd", linestyle="--")
    ax.axvline(geom["tfsf_x"], color="0.4", linestyle=":", label="TFSF")
    ax.axvline(3.0, color="0.2", linestyle="--", label="Debye")
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", nargs="?", type=Path, default=DEFAULT_CASE)
    parser.add_argument(
        "--pvd",
        type=Path,
        default=REPO / "exports/ParaView/single-core/1D_Debye.msh/1D_Debye.msh.pvd",
    )
    parser.add_argument(
        "--png",
        type=Path,
        default=REPO / "exports/SimulationData/single-core/1D_Debye/debye_1d_compare.png",
    )
    args = parser.parse_args()
    case = load_case(args.case)
    geom = build_geom(case)
    print("1D Debye reference, matched to", args.case)
    print(f"  order {case['order']}, elements {len(geom['j'])}, dt {case['dt']}, t_end {case['t_end']}")
    print(f"  eps_inf {case['eps_inf']}, eps_s {case['eps_s']}")
    print(f"  JSON tau {case['tau_json']} s -> solver tau {case['tau']:.6g} light-meters")
    print(f"  Gaussian spread {case['spread']}, mean {geom['mean']:.6g}, TFSF x {geom['tfsf_x']}")
    if case["tau"] > 10.0 * case["t_end"]:
        print("  The pole does not relax during this run: P stays near 0 and the slab acts as eps_inf.")
    times, frames, dt, n_steps = run_reference(case, geom)
    print(f"  integrated {n_steps} steps, stored {len(frames)} frames, dt {dt}")
    t_check, rel = vacuum_check(geom, times, frames, case["spread"])
    print(f"  reference vs analytic right-going pulse at t={t_check:.3f}: relative L2 {rel:.3e}")
    if args.pvd.is_file():
        compare(case, geom, times, frames, args.pvd, args.png)
    else:
        print(f"No ParaView file at {args.pvd}; reference only.")


if __name__ == "__main__":
    main()
