# 3D_RAMC2_50km_Gupta_1GHz

3D DGTD RCS case of the RAM-C II vehicle at 50 km, Mach 9.5 (Eilmer "sub-GHz
targeted" case), with the plasma sheath taken from plasma_kit H1 grid
functions.

## Contents

| File | Role |
|---|---|
| `3D_RAMC2_50km_Gupta_1GHz.msh` | order-2 tetrahedral mesh (MSH 2.2, all elements saved, physical == elementary tags) |
| `3D_RAMC2_50km_Gupta_1GHz.json` | DGTD case: global operator, CUDA run |
| `3D_RAMC2_50km_Gupta_1GHz.rcs.json` | offline `opensemba_rcs` sweep, 0.2-1.5 GHz, backscatter |
| `omega_p.gf` / `nu_c.gf` | plasma_kit DG grid functions on this mesh and FES (solver DoFs) |
| `description.json` | plasma_kit projection metadata (mesh, order, basis, refinement) |

## Mesh layout

Volumes (outside -> inside): 1 = SMA sphere -> NTFF box (vacuum), 2 = NTFF
box -> TFSF box (vacuum), 3 = TFSF box -> Lorentz surface (vacuum, total
field), 4 = flow field (Lorentz).

Surfaces: 9 = PEC (revolved RAM-C wall + base disk, hollow body), 10 = SMA
(**sphere** of radius 1.816 m centred at x = 0.6355 m), 11-16 = TFSF box
faces, 21-26 = NTFF (rcssurface) box faces, 20 = internal Lorentz/vacuum
interface (not referenced by the JSON).

The TFSF and NTFF interfaces are axis-aligned rectangular boxes (planar to
1e-13 m). The SMA absorber is a sphere, following the 3D RCS reference cases
(`3D_RCS_PEC_1m`, `3D_RCS_Sphere_Box_1m_*`): a sphere gives near-normal
incidence everywhere, while an SMA box reflects at grazing incidence on its
sides (about -20 dB, see `docs/designing-a-box-pml.md`).

Sizing: interest frequency = max f_p of the archive = 5.703e9/(2*pi) =
9.077e8 Hz. `h = c/(5 f) = 66 mm` at and inside the TFSF box (order-2
Gauss-Lobatto, 5 elements per wavelength). The plasma layer uses the largest
size that fits its local wall-to-shock thickness, `h = min(66 mm, gap(x))`,
i.e. one element across the layer, isotropic in the other directions. The
absorber relaxes to about `c/(2.3 f) = 130 mm`.

Result: 123,401 tetrahedra, 171,692 nodes. `minSICN` per volume: 1 = 0.34,
2 = 0.30, 3 = 0.14, 4 = 0.23 (medians 0.82, 0.70, 0.80, 0.82); no negative
Jacobians. Smallest edge 13.3 mm, so the default time step is 7.97e-4
light-metres.

## Lorentz material

Volume 4 uses `gridfunction` profiles pointing at `omega_p.gf` and
`nu_c.gf` (SI rad/s; only 1/c0 applied). The files are written by
plasma_kit `mode: "dgtd"` directly in the DoFs of this mesh and
`DG_FECollection(order=2, basis_type=1)`. DGTD loads them once at startup
and distributes them over MPI; it squares `omega_p` at the quadrature
points. A positivity limiter in plasma_kit keeps the projected rates
non-negative. `f_1 = 0` (cold plasma / Drude).

## RAM / run mode

This case is configured for the assembled `global` operator on CUDA
(`runmode: "cuda-1"`). Note that the earlier 143k-tet version of this case
exceeded a 16 GB GPU. This mesh is 123k tets (14% fewer), but the global
operator may still not fit. **Prefer
`3D_RAMC2_50km_Gupta_1GHz_hesthaven`**, which is the same mesh with the
matrix-free operator.

## Run

From the repository root:

```sh
mpiexec -n 1 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i testData/testInputs/3D_RAMC2_50km_Gupta_1GHz/3D_RAMC2_50km_Gupta_1GHz.json
mpiexec -n 1 ./build/gnu-release-cuda/bin/opensemba_rcs \
  -i testData/testInputs/3D_RAMC2_50km_Gupta_1GHz/3D_RAMC2_50km_Gupta_1GHz.rcs.json
```

## Regeneration

```sh
python3 plasma_kit/src/scripts/post_mesher/post_mesher.py \
  --input plasma_kit/src/cfd/ns_solver/cases/ramc2_subghz_targeted/lmrsim/vtk/fluid-s0012.pvtu \
  --frequency 9.07695147e8 --ppw-out 2.3 \
  --plasma-archive plasma_kit/output/ramc2_subghz_targeted_1255_2d_gupta \
  --evolution-operator global --runmode cuda-1 \
  --rcs-start 2.0e8 --rcs-end 1.5e9 --rcs-steps 141 \
  --out-dir testData/testInputs/3D_RAMC2_50km_Gupta_1GHz
```

Then create the plasma grid functions in the solver's DoFs:

```sh
plasma-fields --config /path/to/preprocess.json
```

with a config such as (relative paths resolve from the config file's
directory; adapt `input` and `mesh` to your checkout layout):

```json
{
  "mode": "dgtd",
  "input": "../../plasma_kit/src/cfd/ns_solver/cases/ramc2_subghz_targeted/lmrsim/vtk/fluid-s0012.pvtu",
  "closure": "gupta",
  "dgtd": {
    "mesh": "3D_RAMC2_50km_Gupta_1GHz.msh",
    "case_json": "3D_RAMC2_50km_Gupta_1GHz.json"
  }
}
```

The 2D archive is only used by the post-mesher to estimate the plasma
rates for the time step; it is not shipped in the case.
