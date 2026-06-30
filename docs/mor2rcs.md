# MOR to RCS post-processing

`opensemba_mor2rcs` replays MOR full-state snapshots (`xfull/x_k`) into the on-disk `rcssurface` probe format, then runs the existing `RCSSurfacePostProcessor` unchanged.

## Expected folder layout

Run from the MOR reconstruction case directory:

```text
.
├── <case>.json
├── <mesh>.msh
└── xfull/
    ├── x_0
    ├── x_1
    └── ...
```

Snapshot format matches `mor_state` export and `opensemba_mor2paraview`:

1. First line: absolute simulation time
2. Second line: vector size (`6 × ndofs`)
3. Remaining lines: `[Ex | Ey | Ez | Hx | Hy | Hz]` packed DOFs

## Usage

```sh
./build/gnu-release-mpi/bin/opensemba_mor2rcs \
  --xdir ./xfull \
  --out ./rcs_output \
  -i testData/rcsInputs/2D_RCS_MOR_Single.rcs.json
```

With explicit paths:

```sh
./build/gnu-release-mpi/bin/opensemba_mor2rcs \
  --case ./my_case.json \
  --mesh ./my_mesh.msh \
  --xdir ./xfull \
  --tags 2 3 4 5 \
  --name mor_rcs \
  --out ./rcs_output \
  -i ./rcs_sweep.json
```

Defaults when run in a case folder:

- Auto-detect single `.json` and `.msh` in cwd
- `--xdir`: `./xfull` if present, else `./x`
- `--out`: `./rcs_output`
- `--name`: `mor_rcs`

## RCS sweep JSON (`-i`)

MOR-specific format (no `runmode` / `casename`):

```json
{
  "frequencies": { "start": 1e6, "end": 250e6, "steps": 301 },
  "angles": {
    "theta": { "start": 1.5707963267948966, "end": 1.5707963267948966, "steps": 1 },
    "phi":   { "start": 0.0, "end": 6.28318530718, "steps": 361 }
  },
  "max_time": null
}
```

## NTFF surface tags

Resolved in order:

1. `--tags` on the command line
2. `probes.rcssurface[].tags` in the case JSON
3. `sources[type=planewave].tags`

MOR cases without an `rcssurface` probe typically use planewave source tags.

## Output

```text
rcs_output/
  RCSSurface/mor_rcs/
    rank0/
      mesh
      surface_data.bin
    farfield/farfieldData_Th_<theta>_Phi_<phi>_dgtd.dat
    rcs/rcsData_Th_<theta>_Phi_<phi>_dgtd.dat
```

The `rank0/` subtree matches live `rcssurface` simulation export and can be fed to `opensemba_rcs` for separate frequency/angle sweeps.

## Notes

- Single-rank execution only.
- 2D RCS requires `theta = π/2` in the sweep JSON.
- Plane-wave incident power for normalization is taken from the case JSON `sources` block.
