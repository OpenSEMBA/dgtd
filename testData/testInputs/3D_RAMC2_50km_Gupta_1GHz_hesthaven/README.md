# 3D_RAMC2_50km_Gupta_1GHz_hesthaven

Same 50 km RAM-C II case as `3D_RAMC2_50km_Gupta_1GHz`, but configured for the
matrix-free **hesthaven** operator (CPU/MPI, no assembled global matrix). The
`global` variant exceeded a 16 GB GPU on the earlier mesh, so this is the
variant to run.

The hesthaven operator rejects TFSF faces that touch curved elements
(`HesthavenEvolution::checkForTFSFInCurvedElements`). With the current sizing
the TFSF box needs **no widening**: the post-mesher's curvature check, which
reproduces DGTD's classifier on the written MSH file, reports **0** TFSF
faces touching a curved element, and DGTD's own
`initCurvedAndLinearElementsLists` agrees (0 bad faces). The
`--ensure-straight-tfsf` flag remains available if a future mesh needs it.

## Contents

| File | Role |
|---|---|
| `3D_RAMC2_50km_Gupta_1GHz_hesthaven.msh` | 123,401 order-2 tets, 171,692 nodes |
| `3D_RAMC2_50km_Gupta_1GHz_hesthaven.json` | hesthaven case, gridfunction Lorentz |
| `3D_RAMC2_50km_Gupta_1GHz_hesthaven.rcs.json` | 0.2-1.5 GHz backscatter sweep, runmode single-core |
| `plasma/` | plasma_kit H1 archive (mesh.mesh, omega_p.gf, nu_c.gf) |

Mesh layout, sizing and Lorentz material are identical to the `global` case
(see its README): interest frequency = max f_p = 9.077e8 Hz, 5 elements per
wavelength at/inside the TFSF box, plasma layer at `min(66 mm, gap(x))`,
hollow PEC body (tag 9), TFSF box 11-16, NTFF box 21-26, **spherical** SMA
(tag 10, R = 1.816 m), Lorentz volume 4.

`minSICN` per volume: 1 = 0.34, 2 = 0.30, 3 = 0.14, 4 = 0.23 (medians 0.82,
0.70, 0.80, 0.82); no negative Jacobians. Smallest edge 13.3 mm, default time
step 7.97e-4 light-metres.

## Run

From the repository root:

```sh
./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i testData/testInputs/3D_RAMC2_50km_Gupta_1GHz_hesthaven/3D_RAMC2_50km_Gupta_1GHz_hesthaven.json
./build/gnu-release-mpi/bin/opensemba_rcs \
  -i testData/testInputs/3D_RAMC2_50km_Gupta_1GHz_hesthaven/3D_RAMC2_50km_Gupta_1GHz_hesthaven.rcs.json
```

(or `mpiexec -n N` for MPI). If the machine is CUDA-capable and you prefer
hesthaven on GPU, run the cuda build; `runmode` in the rcs.json would then be
`cuda-N`.

## Regeneration

```sh
python3 plasma_kit/scripts/post_mesher/post_mesher.py \
  --input plasma_kit/ns_solver/cases/ramc2_subghz_targeted/lmrsim/vtk/fluid-s0012.pvtu \
  --frequency 9.07695147e8 --ppw-out 2.3 \
  --plasma-archive testData/testInputs/3D_RAMC2_50km_Gupta_1GHz_hesthaven/plasma \
  --evolution-operator hesthaven --runmode single-core \
  --ensure-straight-tfsf --tfsf-widen-attempts 4 \
  --rcs-start 2.0e8 --rcs-end 1.5e9 --rcs-steps 141 \
  --out-dir testData/testInputs/3D_RAMC2_50km_Gupta_1GHz_hesthaven
```

`--ensure-straight-tfsf` widens the TFSF box by 0.5*h per attempt and remeshes
until no TFSF face touches a curved element (required by hesthaven). Use
`--mesh-order 1` instead if you want a fully straight geometry (G1_O2).
