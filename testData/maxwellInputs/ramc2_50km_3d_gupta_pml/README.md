# RAM-C II 50 km, 3D Gupta, PML absorber

Prepared case for DGTD.  The plasma rates are DG grid functions produced by
plasma_kit `mode.solver: "OpenSEMBA-dgtd"` on this mesh and finite element space (order-2 H1
source projected onto the order-2 DG space), stored beside the case JSON.

## Contents

| File | Role |
|---|---|
| `ramc2_50km_3d_gupta_pml.msh` | 3D order-2 tetrahedral mesh, 118,678 tets, 169,656 nodes |
| `ramc2_50km_3d_gupta_pml.json` | DGTD case, `global` operator, PML absorber |
| `ramc2_50km_3d_gupta_pml.rcs.json` | offline RCS sweep, 0.2–1.5 GHz, backscatter, `runmode: "mpi-16"` |
| `ramc2_50km_3d_gupta_pml.rcs.cuda-1.json` | same sweep with `runmode: "cuda-1"` for single-rank GPU runs |
| `omega_p.gf` / `nu_c.gf` | plasma_kit DG grid functions in the solver's DoFs |
| `description.json` | plasma_kit projection metadata (mesh, order, basis) |

## Material and tag layout

Outside to inside:

| Tags | Region | Material |
|---|---|---|
| 1..26 | Cartesian PML layer (X/Y/Z faces, XY/XZ/YZ edges, XYZ corners) | `PML`, `grading_order: 2` |
| 27 | NTFF shell (NTFF box -> TFSF box) | vacuum |
| 28 | TFSF shell (TFSF box -> flow field) | vacuum |
| 29 | flow field | `lorentz` with `gridfunction` `f_p` / `gamma` files |

Boundaries: `9` = PEC, `10` = SMA on the PML outer truncation, `11..16` =
TFSF faces, `21..26` = NTFF faces.  `20` is the internal Lorentz/vacuum
interface and is not referenced by the JSON.

PML blocks use `grading_order: 2`, `target_reflection: 1e-6`,
`kappa_max: 1.0`, and `alpha_max: 0.0`.  Face tags are 1..6, edge tags 7..18,
and corner tags 19..26.

## Run

From the directory above the case:

```sh
opensemba_dgtd -i ramc2_50km_3d_gupta_pml/ramc2_50km_3d_gupta_pml.json
opensemba_rcs  -i ramc2_50km_3d_gupta_pml/ramc2_50km_3d_gupta_pml.rcs.json
```

The default `.rcs.json` is set to `runmode: "mpi-16"`.  For a single-rank
CUDA run, use `ramc2_50km_3d_gupta_pml.rcs.cuda-1.json` instead.

Key run values:

- DG order 2, Gauss-Lobatto basis,
- `global` operator, TFSF curvature check passed (0/6716 bad faces),
- time step `8.17491e-4` light-metres,
- final time `23.6692` light-metres,
- `rcssurface` sampling every 13 steps,
- exporter `context_fields: "all"`, `context_mode: "auto"` (first frame only).

## Transfer

From this machine, for example:

```sh
scp -r output/ramc2_50km_3d_gupta_pml user@host:/path/to/plasma_kit/output/
```

The case is self-contained once the folder is copied; the JSON points at
`omega_p.gf` and `nu_c.gf` beside it.
