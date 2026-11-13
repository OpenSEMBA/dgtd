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
| `plasma/` | plasma_kit archive: `mesh.mesh` + `omega_p.gf` + `nu_c.gf` + `description.json` |

## Mesh layout

Volumes (outside -> inside): 1 = SMA->NTFF (vacuum), 2 = NTFF->TFSF
(vacuum), 3 = TFSF->Lorentz surface (vacuum, total field), 4 = flow field
(Lorentz). Surfaces: 9 = PEC (revolved RAM-C wall + base disk, hollow body),
10 = SMA, 11-16 = TFSF box faces, 21-26 = NTFF (rcssurface) faces, 20 =
internal Lorentz/vacuum interface (not referenced by the JSON).

Sizing: interest frequency = max f_p of the archive = 5.703e9/(2*pi) =
9.077e8 Hz; h = c/(5 f) = 66 mm at and inside the TFSF box (order-2
Gauss-Lobatto, 5 elements per wavelength), relaxing outwards, capped at
c/(3 f); 3 elements across the 24 mm plasma layer. 143,142 tetrahedra,
265,710 nodes.

## Lorentz material

Volume 4 uses `gridfunction` profiles pointing at `plasma/omega_p.gf` and
`plasma/nu_c.gf` (SI rad/s; only 1/c0 applied; 2D archive evaluated
axisymmetrically as (x, r=hypot(y,z)); `fallback: "nearest"`). DGTD evaluates
the H1 functions at the DG quadrature points and applies the per-element
M_K^-1 M_K[f] projection. `f_1 = 0` (cold plasma / Drude).

## RAM / run mode

Estimated global-operator RAM ~ 8 GB (< 16 GB), so the case is configured
for `evolution_operator: "global"` on CUDA (`runmode: "cuda-1"` in the
rcs.json). If the target machine reports more than 16 GB, switch to
`"hesthaven"` and `"runmode": "single-core"`.

## Run

From the repository root:

```sh
mpiexec -n 1 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i testData/testInputs/3D_RAMC2_50km_Gupta_1GHz/3D_RAMC2_50km_Gupta_1GHz.json
mpiexec -n 1 ./build/gnu-release-cuda/bin/opensemba_rcs \
  -i testData/testInputs/3D_RAMC2_50km_Gupta_1GHz/3D_RAMC2_50km_Gupta_1GHz.rcs.json
```

## Regeneration

The mesh and JSONs were produced with plasma_kit's post-mesher:

```sh
python3 plasma_kit/scripts/post_mesher/post_mesher.py \
  --input plasma_kit/ns_solver/cases/ramc2_subghz_targeted/lmrsim/vtk/fluid-s0012.pvtu \
  --frequency 9.07695147e8 \
  --plasma-archive testData/testInputs/3D_RAMC2_50km_Gupta_1GHz/plasma \
  --evolution-operator global --runmode cuda-1 \
  --rcs-start 2.0e8 --rcs-end 1.5e9 --rcs-steps 141 \
  --out-dir testData/testInputs/3D_RAMC2_50km_Gupta_1GHz
```

(plasma_kit and dgtd checkouts side by side). The plasma archive originates
from `plasma_kit/output/ramc2_subghz_targeted_1255_2d_gupta`
(Gupta closure, Eilmer 50 km targeted case).
