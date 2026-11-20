# 3D_RAMC2_50km_Gupta_1GHz_hesthaven

Same 50 km RAM-C II case as `3D_RAMC2_50km_Gupta_1GHz`, but configured for the
matrix-free **hesthaven** operator (CPU/MPI, no assembled global matrix). The
`global` variant exceeded a 16 GB GPU on the earlier mesh, so this is the
variant to run.

The hesthaven operator rejects TFSF faces that touch curved elements
(`HesthavenEvolution::checkForTFSFInCurvedElements`). With the current sizing
the TFSF box needs **no widening**: the mesh builder's curvature check, which
reproduces DGTD's classifier on the written MSH file, reports **0** TFSF
faces touching a curved element, and DGTD's own
`initCurvedAndLinearElementsLists` agrees (0 bad faces). The
`meshing.ensure_straight_tfsf` option remains available if a future mesh needs
it.

## Contents

| File | Role |
|---|---|
| `3D_RAMC2_50km_Gupta_1GHz_hesthaven.msh` | 123,401 order-2 tets, 171,692 nodes |
| `3D_RAMC2_50km_Gupta_1GHz_hesthaven.json` | hesthaven case, gridfunction Lorentz |
| `3D_RAMC2_50km_Gupta_1GHz_hesthaven.rcs.json` | 0.2-1.5 GHz backscatter sweep, runmode single-core |
| `omega_p.gf` / `nu_c.gf` | plasma-kit DG grid functions on this mesh and FES |
| `description.json` | plasma-kit projection metadata |

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

Create `plasmakit-generate.json` (for example in the repository root) with:

```json
{
  "mode": {
    "solver": "OpenSEMBA-dgtd",
    "build": {
      "output_mesh": "testData/testInputs/3D_RAMC2_50km_Gupta_1GHz_hesthaven/3D_RAMC2_50km_Gupta_1GHz_hesthaven.msh",
      "frequency": 9.07695147e8,
      "order": 2,
      "basis_type": 1,
      "overwrite": true,
      "sizing": {"ppw_out": 2.3},
      "boxes": {"absorber": "sma"},
      "physics": {
        "evolution_operator": "hesthaven",
        "timestep_archive": "plasma_kit/output/ramc2_subghz_targeted_1255_2d_gupta"
      },
      "rcs": {"start": 2.0e8, "end": 1.5e9, "steps": 141,
              "runmode": "single-core"},
      "meshing": {"ensure_straight_tfsf": true, "tfsf_widen_attempts": 4}
    }
  },
  "input": "plasma_kit/src/cfd/ns_solver/cases/ramc2_subghz_targeted/lmrsim/vtk/fluid-s0012.pvtu",
  "plasma_model": "gupta"
}
```

Then run:

```sh
plasma_kit/build/gnu-release/bin/plasmakit --config plasmakit-generate.json
```

The hesthaven mesh auto-enables TFSF widening when `mesh_order > 1`; the
explicit `ensure_straight_tfsf: true` is kept here for clarity. Use
`mesh_order: 1` instead if you want a fully straight geometry (G1_O2).

`timestep_archive` is only a time-step reference (it supplies
`omega_p_max` / `nu_c_max`); the fields are recomputed from `input`.
