# Designing a box PML

Cartesian SC-PML (Bagci/Chen) on `GlobalEvolution`. The mesh owns the layer geometry. JSON only names which element tags are PML and how conductivity is graded inside them. There is no thickness field.

Formulation: [30-sc-pml-ade.md](./30-sc-pml-ade.md). The absorption numbers below are the reference for this absorber.

Units inside the solver are normalized: $c=\varepsilon=\mu=1$. Mesh coordinates are those units. Use `"evolution_operator": "global"` or omit it. `"hesthaven"` does not run volumetric PML.

## Cut the box

Draw the domain of interest, then a layer of PML around it, then stop the mesh. The vacuum–PML faces must be planes perpendicular to the coordinate axes, and each side must sit at one coordinate. The solver finds those planes from faces where a vacuum element touches a PML element, and grades depth from that plane. A slanted, wavy, or stair-stepped cut is graded from a single plane anyway, so the conductivity will not follow the elements you drew.

Split the layer so each element has a fixed set of stretch axes:

| Region | Where | `active_axes` |
|--------|--------|----------------|
| Interior | Everything inside the cut | not PML (`"type": "vacuum"`) |
| Face slab | Touches one outer face only. In 2D that is a left/right or top/bottom strip, corners excluded | the outward axis only: `"X"`, `"Y"`, or `"Z"` |
| Edge | Two slabs meet. In 2D these are the four corners | the two axes, e.g. `["X","Y"]` |
| Corner | Three slabs meet (3D only) | `["X","Y","Z"]` |

A 2D box is four slabs plus four corners. A 3D box is six face slabs, twelve edges, and eight corners. Opposite slabs that share an axis can share one JSON block (two tags). Corners and edges need their own block: an element tagged as an $X$-only slab is not stretched in $Y$ or $Z$.

The gold-standard mesh does this in 2D. Vacuum occupies $|x|,|y| \lesssim 3$ (tags 9 and 10). PML fills out to the outer box at about 4.

| Tags | Role | `active_axes` |
|------|------|----------------|
| 1, 2 | top and bottom slabs | `["Y"]` |
| 3, 4 | left and right slabs | `["X"]` |
| 5, 6, 7, 8 | four corners | `["X","Y"]` |
| 9, 10 | interior | vacuum |

In Gmsh, material tags are physical groups of the elements (surfaces in 2D, volumes in 3D). Boundary tags are separate physical groups of the boundary entities (curves in 2D, surfaces in 3D). Put **SMA** on the outer boundary of the PML. The PML is the volume; SMA is only the truncation behind it. Do not use PEC there for an absorption test.

Each element attribute may appear in only one PML block. A tag that is not in the mesh, or that is listed twice, is an error at startup.

No 3D box PML case is in the tree yet. The same cut applies: one plane per side, face / edge / corner tags as in the table, SMA on the outer surface. Judge a new 3D mesh with the measurement in the last section. Do not expect the 2D decibel numbers on a different geometry.

## What the solver does with the cut

For each axis that has a vacuum–PML face, depth at a point $x$ is how far that coordinate sits past the interface plane. On the $+X$ side that is $x - x_{+}$. On the $-X$ side it is $x_{-} - x$. The same rule applies to $Y$ and $Z$.

Thickness $L$ is not in the JSON. For each PML block and each of its active axes, $L$ is the largest such depth among quadrature points in that block. On the dipole mesh, $L \approx 1$.

Normalized depth is $\xi = \rho / L$, clamped to $[0,1]$. With grading order $m$ and design reflection $R$ (`target_reflection`):

$$
\sigma_{\max} = -\frac{(m+1)\ln R}{2L}
$$

$$
\sigma(\xi) = \sigma_{\max}\,\xi^{m}
$$

$m = 0$ means $\xi^{0} = 1$, so $\sigma$ is the constant $\sigma_{\max}$ through the layer. $R$ sets that constant. It is a single-interface design estimate. The measured probe reflection in the table below is weaker than $20\log_{10}(R)$. For $R = 10^{-6}$ that estimate is $-120\,\mathrm{dB}$; the dipole measurement is about $-80\,\mathrm{dB}$.

Inactive axes in a block stay at $\sigma = 0$. A corner block grades each of its axes from its own plane, using that block's own $L$.

## JSON attributes

One object per tag group, inside `model.materials`. `active_axes` is required. Everything else has a default.

| Field | Default | What you can set |
|-------|---------|------------------|
| `tags` | — | Element attribute ids for this block. Required. |
| `type` | — | `"PML"`. |
| `active_axes` | — | One or more of `"X"`, `"Y"`, `"Z"` (case-insensitive). At least one. An axis past the mesh dimension is an error. `"R"` is an error. |
| `matches_vacuum` | `true` | Must be `true` or omitted. The layer is free space plus stretch. `false` is an error. |
| `grading_order` | `3` | Integer $m \ge 0$. `0` is constant $\sigma$. Larger $m$ holds $\sigma$ down near the interface and piles it up at the outer boundary. On the reference mesh, `0` absorbed best on the axis probes. |
| `target_reflection` | `1e-6` | Design $R$ in $(0,1)$. Smaller $R$ raises $\sigma_{\max}$. This is not the decibel number you will measure. |
| `kappa_max` | `1` | $\kappa \ge 1$ at the outer boundary. At $1$, $\kappa \equiv 1$ everywhere (the reference runs). Above $1$, $\kappa(\xi) = 1 + (\kappa_{\max}-1)\,\xi^{m}$, and $m=0$ uses $\kappa_{\max}$ as a constant. |
| `alpha_max` | `0` | Must be `0` or omitted. A CFS pole is not implemented. |
| `stretch_mode` | `"box"` | `"box"` or `0`. Depth is the planar depth above. `"radial"` / `1` still parses and grades by radius, but radial PML cases have been removed. New meshes use `"box"`. |
| `radial_center` | — | Only legal with `stretch_mode` `"radial"`. Do not set it on a box. |

Rejected on a PML block: `sigma_max`, `bulk_conductivity`, `relative_permittivity`, `relative_permeability`. Conductivity of the layer comes only from $R$, $m$, and $L$.

The reference case `testData/maxwellInputs/2D_Dipole_PML_GO0/2D_Dipole_PML_GO0.json` is the block to copy. `GO1`, `GO2`, and `GO3` change only `grading_order`. The SMA twin `2D_Dipole` uses the same outer box with every tag set to vacuum.

```json
{
  "tags": [3, 4],
  "type": "PML",
  "matches_vacuum": true,
  "grading_order": 0,
  "target_reflection": 1e-6,
  "kappa_max": 1.0,
  "alpha_max": 0.0,
  "active_axes": ["X"]
}
```

## Absorption reference

This is the measurement a box PML is compared against. Same source, order 3, $\Delta t = 0.005$, $T = 16$, upwind 1, MPI 4. Probes sit in the vacuum, inside the interface: X $(2.5, 0)$, Y $(0, 2.5)$, XY $(1.768, 1.768)$. Component $E_z$. Incident window about $[3.27, 4.08]$, reflection window about $[6.46, 16]$ (`--ref-delay 3`).

$$
R_{\mathrm{dB}} = 20\log_{10}\big(|E_{\mathrm{ref}}| / |E_{\mathrm{inc}}|\big)
$$

at the incident spectral peak. Numbers below are the reproduced run (they match [31-sc-pml-status.md](./31-sc-pml-status.md) §8 when rounded).

| Case | $m$ | X | Y | XY |
|------|-----|----|----|-----|
| `2D_Dipole` (SMA only) | — | −20.00 | −20.00 | −30.53 |
| `2D_Dipole_PML_GO0` | 0 | −79.69 | −82.58 | −76.81 |
| `2D_Dipole_PML_GO1` | 1 | −72.39 | −72.59 | −82.29 |
| `2D_Dipole_PML_GO2` | 2 | −74.40 | −72.12 | −70.15 |
| `2D_Dipole_PML_GO3` | 3 | −74.03 | −70.88 | −73.69 |

SMA alone sits near $-20\,\mathrm{dB}$ on axis. The PML mesh, same box, drops that by about $50$–$60\,\mathrm{dB}$. Flat grading ($m = 0$) is the best on-axis result. Higher $m$ is a few decibels worse on X/Y and, for $m = 1$, better on the diagonal.

Replay after a run. Point probes append, so remove the case export directory before repeating a case. The `GO*` JSON files share one mesh file; run those cases one after another, not in one `mpirun`.

```sh
python3 scripts/dipole_pml_grading_dft.py \
  --ref-delay 3 \
  Exports/SimulationData/mpi-4/2D_Dipole \
  Exports/SimulationData/mpi-4/2D_Dipole_PML_GO0 \
  Exports/SimulationData/mpi-4/2D_Dipole_PML_GO1 \
  Exports/SimulationData/mpi-4/2D_Dipole_PML_GO2 \
  Exports/SimulationData/mpi-4/2D_Dipole_PML_GO3
```

A new box is the same kind of comparison: an all-vacuum SMA twin of that outer boundary, probes in the vacuum, and this DFT. Quote $R_{\mathrm{dB}}$ next to the table. Matching $-80\,\mathrm{dB}$ is what this mesh and this grading achieved, not a guarantee for a thicker layer, a coarser mesh, or a 3D box.
