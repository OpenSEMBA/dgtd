# 2D hybrid structured / unstructured box + circle

Demonstration mesh for Hesthaven-friendly layout: **structured triangles** (transfinite lattice) in the TFSF / NTFF / PML frame, **unstructured triangles** only between the circle and the TFSF box. Transfinite corners are CCW so element Jacobians stay positive; `remesh.sh` also flips any leftover inverted tris.

## Geometry

| Region | Half-size | Role |
|--------|-----------|------|
| Circle | $r=0.5$ | Object / disk (surface tag 10, curve tag 1) |
| TFSF box | $1$ | Closed TF/SF interface (curve tag 9) |
| NTFF box | $3$ | Huygens / `rcssurface` (curve tag 11) |
| PML interface | $4$ | Vacuum–PML cut |
| Outer box | $5$ | SMA truncation |

Domain: $[-5,5]^2$.

## Meshing idea

- **Core** (`|x|,|y|\le 1` minus disk): Frontal-Delaunay triangles.
- **Rings** $1\to 3$, $3\to 4$, $4\to 5$: `Transfinite` **triangles** (structured lattice, two tris per logical quad). No `Recombine` — MFEM/`opensemba_dgtd` requires a single element type in the mesh.
- Every concentric square uses the **same tangential node count** (`n_tang = 80`). That is required so opposite sides of each transfinite patch match. Radial counts (`n_r13`, `n_r34`, `n_r45`) refine only in the radial direction.

### Outer → inner meshing (faces)

Topology is built inside-out (shared edges), but **surfaces are meshed outer → inner**:

1. PML ($4\to 5$)
2. Vacuum ring $3\to 4$
3. Vacuum ring $1\to 3$
4. TFSF core
5. Disk

Elementary surface tags are numbered that way too (1…4 = PML, …, 14 = disk). This keeps shared-face / Elem1–Elem2 orientation anchored from the truncated exterior. Do not remesh with bare `gmsh -2` (that can reorder); use `./remesh.sh`.

## Surface tags (materials)

| Tag | Region |
|-----|--------|
| 9 | Vacuum (core + structured rings inside PML) |
| 10 | Disk |
| 1 | PML $+Y$ slab |
| 2 | PML $-Y$ slab |
| 3 | PML $-X$ slab |
| 4 | PML $+X$ slab |

**Note:** PML corners are folded into the four face quads (one cell deep). Strict Dipole_PML corner tags 5–8 would need an 8-block outer ring (mid-side points on $|\,|=4$ and $|\,|=5$). Face-only `active_axes` is still valid for a first hybrid demo.

## Boundary tags

| Tag | Role |
|-----|------|
| 1 | Circle |
| 9 | TFSF |
| 11 | NTFF / RCS surface |
| 4, 6, 13, 15 | Outer SMA |

## Size (default `.geo`)

Rebuild prints ~30k nodes / ~35k elements (order-1). Raise `n_tang` / `n_r*` in the `.geo` for a heavier mesh.

## Rebuild

```sh
./remesh.sh
```

Requires Gmsh ≥ 4.x on `PATH`. The optional `build_mesh.py` is a sketch of the same topology via the Gmsh Python API (not required if the `.geo` builds).

## Sample JSON

`2D_Hybrid_Box_Circle.json` — Hesthaven + planewave on TFSF tag 9, PML face slabs, NTFF `rcssurface` tag 11. Tune `time_step` / `order` before production runs.
