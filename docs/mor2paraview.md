# MOR to ParaView post-processing

`opensemba_mor2paraview` replays previously exported `mor_state` snapshots (`x_0`, `x_1`, …) into a ParaView time series dataset.

## Expected folder layout

Run the command from the case directory:

```text
.
├── <case>.json
├── <mesh>.msh
└── x/
    ├── x_0
    ├── x_1
    └── ...
```

## Snapshot format (`x_k`)

1. First line: absolute simulation time  
2. Second line: vector size  
3. Remaining lines: state values in `[Ex | Ey | Ez | Hx | Hy | Hz]` packed order  

## Usage

Default (auto-detect single `.json` and `.msh` in cwd):

```sh
./build/gnu-release-mpi/bin/opensemba_mor2paraview
```

With explicit paths:

```sh
./build/gnu-release-mpi/bin/opensemba_mor2paraview \
  --case ./my_case.json \
  --mesh ./my_mesh.msh \
  --xdir ./x \
  --out ./output \
  --name mor_paraview.vtk
```

Defaults:

- Reads snapshots from `./x`
- Writes to `./output`
- Collection name `mor_paraview.vtk`

## Notes

- Single-rank execution only.
- Output format matches the exporter probe time-series layout for ParaView.
