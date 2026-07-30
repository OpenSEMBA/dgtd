# MOR to RCS post-processing

`opensemba_mor2rcs` builds the on-disk `rcssurface` probe format and runs
`RCSSurfacePostProcessor`. Two snapshot sources are supported:

1. **Legacy `--xdir`** — ASCII full-order MOR dumps (`x_k`, size `6×ndofs`).
2. **Reconstruct `--ur` + `--xrdir`** — dense basis `Ur.bin` and reduced states
   `xr_k`; full-order state is formed on the fly as `x = Ur * xr` and **never**
   written unless `--dump-xdir` is set.

Peak memory in reconstruct mode is about **`Ur` + one full `x`** (almond-scale:
~tens of GiB for `Ur`, ~45 MiB for `x`). A multi-terabyte ASCII `xfull/` tree is
**not** required.

## Expected folder layouts

### Legacy (`--xdir`)

```text
.
├── <case>.json
├── <mesh>.msh
└── xfull/
    ├── x_0
    ├── x_1
    └── ...
```

Snapshot format (same as `mor_state` / `mor2paraview`):

1. First line: absolute simulation time  
2. Second line: vector size (`6 × ndofs`)  
3. Remaining lines: `[Ex | Ey | Ez | Hx | Hy | Hz]` packed DOFs  

### Reconstruct (`--ur` / `--xrdir`)

Produced by MOR (`mor_dgtd`) when operators are saved:

```text
Exports/<mor_case>/
  operators/
    Ur.bin              # dense float64, column-major (N × r), host endian
    meta.json           # optional but recommended
    Ar.bin / Br.bin     # ignored by mor2rcs
  xr/
    xr_0.txt            # or bare xr_0; time, size(=r), then r coeffs (ASCII)
    xr_1.txt
    ...
```

`xr` files may be named `xr_<k>` or `xr_<k>.txt` (MOR uses the `.txt` suffix).

#### `meta.json` contract

```json
{
  "N": 5619180,
  "r": 282,
  "dtype": "float64",
  "layout": "colmajor",
  "ur_file": "Ur.bin",
  "dt": 0.0005,
  "final_time": 16.0,
  "mor_case": "3D_Nasa_Almond_G2_25cm_5GHz_NEWS_train_Base_evol"
}
```

- **`layout`**: only `colmajor` (Fortran): entry `(i,j)` at byte offset
  `(i + j*N) * 8`.  
- **`dtype`**: only `float64` / `double`.  
- If `meta.json` sits next to `Ur.bin`, it is loaded automatically; CLI `--meta`
  overrides. Without meta, `N` is taken from the FE space and `r` from the first
  `xr_*` header.  
- `Ur` rows **must** equal `fields.fieldBlockSize()`; `Ur` cols **must** equal
  each `xr` size.

## Usage

Legacy:

```sh
./build/gnu-release-mpi/bin/opensemba_mor2rcs \
  --xdir ./xfull \
  --out ./rcs_output \
  -i testData/rcsInputs/2D_RCS_MOR_Single.rcs.json
```

Reconstruct-on-the-fly (no `xfull`):

```sh
./build/gnu-release-mpi/bin/opensemba_mor2rcs \
  --case ./my_case.json \
  --mesh ./my_mesh.msh \
  --ur ./operators/Ur.bin \
  --xrdir ./xr \
  --tags 11 12 13 14 15 16 \
  --name almond \
  --out ./rcs_output \
  -i ./sweep.rcs.json
```

Optional: `--meta ./operators/meta.json`, `--dump-xdir ./x_debug` (write
reconstructed ASCII `x_k` for comparison). If both `--ur/--xrdir` and `--xdir`
are passed, `--xdir` is ignored (with a warning).

Defaults when run in a case folder:

- Auto-detect single `.json` and `.msh` in cwd  
- `--xdir`: `./xfull` if present, else `./x` (legacy mode only)  
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

The `rank0/` subtree matches live `rcssurface` simulation export and can be fed
to `opensemba_rcs` for separate frequency/angle sweeps.

## Notes

- Single-rank execution only.  
- Reconstruct mode keeps the same Galerkin `x = Ur*xr`, same surface tags, and
  the same `RCSSurfacePostProcessor` as legacy `--xdir`.  
- For bit-comparable checks: run both modes on the same case/tags/mesh and
  compare `surface_data.bin` (or `--dump-xdir` vs `xfull`) to ~1e-10 relative.  
- 2D RCS requires `theta = π/2` in the sweep JSON.  
- Plane-wave incident power for normalization is taken from the case JSON
  `sources` block.  
