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

Volume 4 uses `gridfunction` profiles pointing at `plasma/omega_p.gf` and
`plasma/nu_c.gf` (SI rad/s; only 1/c0 applied; 2D archive evaluated
axisymmetrically as (x, r=hypot(y,z)); `fallback: "nearest"`). DGTD evaluates
the H1 functions at the DG quadrature points and applies the per-element
M_K^-1 M_K[f] projection. `f_1 = 0` (cold plasma / Drude).

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
python3 plasma_kit/scripts/post_mesher/post_mesher.py \
  --input plasma_kit/ns_solver/cases/ramc2_subghz_targeted/lmrsim/vtk/fluid-s0012.pvtu \
  --frequency 9.07695147e8 --ppw-out 2.3 \
  --plasma-archive testData/testInputs/3D_RAMC2_50km_Gupta_1GHz/plasma \
  --evolution-operator global --runmode cuda-1 \
  --rcs-start 2.0e8 --rcs-end 1.5e9 --rcs-steps 141 \
  --out-dir testData/testInputs/3D_RAMC2_50km_Gupta_1GHz
```

The plasma archive originates from
`plasma_kit/output/ramc2_subghz_targeted_1255_2d_gupta` (Gupta closure, Eilmer
50 km targeted case).
