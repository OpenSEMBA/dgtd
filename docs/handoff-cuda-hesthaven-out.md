# Handoff: Fix CUDA Hesthaven `out` Vector

**Status:** In progress — sync/layout/init fixes landed; **`out` still zero on GPU** after synced read.  
**Plan file (do not edit):** `.cursor/plans/fix_hesthaven_cuda_out_467a00f6.plan.md`  
**Prior chat:** [Fix CUDA Hesthaven out](10650082-d3af-419c-9598-ebac80812d46)

---

## Goal

CUDA Hesthaven must produce a **nonzero rate vector `out`** so RK4 evolves the field and ParaView shows a visible wave.

**Primary test case** (use exactly this JSON/mesh — do not create alternate case files):

```bash
./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_RCS_PEC_Circle_G1_1m_Uniform_Hesthaven/2D_RCS_PEC_Circle_G1_1m_Uniform_Hesthaven.json \
  --device cuda
```

**Scope:** Hesthaven files only (`HesthavenEvolution.{cpp,cu,h}`). GlobalEvolution is reference-only. TFSF stays **CPU-evaluated** on host jump buffers + `Write()` to device (user choice).

---

## Symptom

| Path | TFSF/jumps | `out` after sync | ParaView E |
|------|------------|------------------|------------|
| CPU (`--device cpu` or no GPU mult) | grows | **nonzero** (~1e-7 early) | visible wave |
| CUDA | grows (host TFSF debug) | **exactly zero** (`stale` and synced) | ~zero |

So: **jumps look alive on host; GPU element pass produces nothing.**

---

## What we've done so far

### 1. `out` sync ordering (`HesthavenEvolution.cpp` — `MultGPU`)

- `out.UseDevice(true); out = 0.0` before kernel
- After `hesthaven_mult_gpu` + curved `AddMult`: `(void)out.Read()` **before** host debug sums
- Final `(void)out.ReadWrite()` for RK (mirrors `GlobalEvolution`)
- Temporary debug: `stale(out)` vs synced `abs(out)` — **both are zero** → not just a stale-host bug

### 2. Static GPU buffer upload (`initGPUData` / `initGPUBoundaryData`)

- Added `hesthaven_sync_gpu_static_data()` — calls `.Write()` on all static device arrays after host fill
- Refactored init to use `HostWrite()` for: matrix metadata, elem ids/dofs, normals/fscale, BC/TFSF index arrays

### 3. Column-major `dense_matvec` (`HesthavenEvolution.cu`)

- Fixed from row-major to Eigen column-major: `mat[i + j * rows]`
- CPU uses `Eigen::MatrixXd` matvec; GPU must match

### 4. Element accumulation path (`HesthavenEvolution.cu`)

- Two-kernel **elem_out + scatter** caused **segfault** → reverted to **direct `out_d[...] +=`** in element kernel (simpler, matches pre-refactor intent)

### 5. Jumps/BC/TFSF path (`MultGPU`)

- Replaced `hesthaven_compute_jumps_gpu` + `hesthaven_apply_bc_gpu` with **full CPU path** (same as `MultCPU`):
  - `computeJumps` → `applyBoundaryConditionsToNodes` → `evaluateTFSF`
  - `memcpy` into `d_jumps_e/h`, then `.Write()`
- **Still zero `out`** → bug is in **`hesthaven_mult_gpu`** or its inputs (`in` on device, metadata, LIFT), not the jump kernels

### 6. Tests run (CUDA build `gnu-release-cuda`)

```bash
./build/gnu-release-cuda/bin/maxwell_solver_tests --device cuda \
  --gtest_filter='HesthavenOperatorBenchmark.GPU_Mult*'
```

| Test | Result | Note |
|------|--------|------|
| `GPU_Mult_matches_CPU_Mult` | PASS | `norm_cpu=0`, `norm_gpu=0` — **trivial match** |
| `GPU_Mult_TFSF_matches_CPU_Mult` | PASS | `||y_cpu - y_gpu||=0` — both may be zero |

Gtests pass but may not exercise nonzero GPU mult (see “Next steps”).

### 7. Out of scope (already done elsewhere)

- `Solver.cpp`: `fields_.allDOFs().HostRead()` before probes on CUDA (ParaView export sync) — keep if still needed

---

## Files touched

| File | Changes |
|------|---------|
| `src/evolution/HesthavenEvolution.cpp` | `hesthaven_sync_gpu_static_data`, init `HostWrite`, `MultGPU` sync + CPU jumps path + debug |
| `src/evolution/HesthavenEvolution.cu` | column-major matvec, direct `out` accumulate, non-const `gpu` in `hesthaven_mult_gpu` |
| `src/evolution/HesthavenEvolution.h` | `hesthaven_mult_gpu(HesthavenGPUData& ...)` signature |

**Not yet done:** remove temporary debug prints (`cleanup-debug` todo).

---

## Leading hypotheses (ordered)

### A. `in` not valid on device when kernel runs (high priority)

- `MultGPU` starts with `(void)in.Read()` (device→host if device newer)
- Jumps use `FieldsInputMaps` → `in.HostRead()` (host path)
- `hesthaven_mult_gpu` uses `in.Read()` (device pointer)
- **Mismatch:** if ODE state or test vectors are host-filled without `UseDevice(true)` + `Write()`, device `in` can be **zero** while host jumps are correct
- **However:** LIFT/flux terms depend on **jumps only**, not `in` — so zero `in` alone may not explain **fully zero** `out` if jumps are correct on device

**Try next:**

```cpp
// After jump host fill (which calls in.HostRead()), before hesthaven_mult_gpu:
in.Write();  // push host→device so kernel sees same state as CPU path
```

Also verify benchmark test sets `x.UseDevice(true)` + `HostWrite()` before `Mult`.

### B. Device jump buffer wrong despite host copy (medium)

- Add diagnostic after `d_jumps_e.Write()`:
  - sum host buffer vs `(void)d_jumps_e.Read();` device readback
- Confirms TFSF actually reached device before element kernel

### C. Element kernel inputs still wrong on device (medium)

- `ref_lift`, `d_matrices`, `d_normals`, `d_fscale`, `d_elem_dofs`, `jump_base = elem_id * flux_size`
- CPU path works with same indexing → if device metadata is zero/stale, kernel silently produces zero
- `hesthaven_sync_gpu_static_data` should have fixed this; **verify** with one-time init checksum print (host vs device sum of `d_matrices`)

### D. Non-atomic `out_d[gi] +=` races (lower)

- Shared face DOFs get concurrent `+=` from multiple elements
- Usually causes **wrong** sums, not **exactly zero everywhere**
- Consider `atomicAdd` only if A–C ruled out

### E. `gpu_.initialized` / `linearElements_.Size()` (ruled out for uniform case)

- Runtime shows **292 linear elements** and element timing — kernel is being called

---

## Plan todos — status

| ID | Task | Status |
|----|------|--------|
| `diag-host-device` | Host vs device diagnostics in `MultGPU` | **Done** (stale vs synced; elem_out vars unused after scatter revert) |
| `fix-out-sync` | `out.Read()` before debug; `ReadWrite()` at end | **Done** |
| `fix-gpu-init-upload` | `.Write()` after init fills | **Done** |
| `fix-matvec-layout` | Column-major `dense_matvec` | **Done** |
| `validate-gtests` | `GPU_Mult_*` on CUDA | **Done** (pass; norms suspicious) |
| `validate-uniform-case` | Uniform circle CUDA → nonzero `out` + ParaView | **Blocked** — still zero |
| `cleanup-debug` | Remove temporary diagnostics | **Pending** (do after fix verified) |

---

## Where we want to get

1. **Device-synced `out` debug grows** on uniform Hesthaven CUDA case (same order of magnitude as CPU early timesteps)
2. **ParaView** shows nonzero E field on CUDA
3. **`GPU_Mult_matches_CPU_Mult`** and **`GPU_Mult_TFSF_matches_CPU_Mult`** pass with **nonzero** `norm_cpu` and small relative error
4. Remove temporary `stale(out)` / `elem_out` debug lines; keep normal timing if desired

---

## Recommended next steps (resume here)

1. **`in.Write()` before `hesthaven_mult_gpu`** in `MultGPU` (after CPU jump path that calls `HostRead`)
2. **Device jump sum diagnostic** — one line comparing host vs device jump norms after `Write()`
3. **One-shot Mult compare** on uniform case: force `benchmarkMultCPU` vs `Mult` with `x.UseDevice(true)`, print `||y_cpu - y_gpu||` and norms
4. If still zero: **single-element CPU vs GPU** — pick `e=0`, print elem flux / LIFT / scatter contributions on CPU; add matching device readback for same elem (or temporary kernel printf)
5. Fix root cause in `HesthavenEvolution.cu` / `.cpp` only
6. Re-run uniform case + gtests
7. **`cleanup-debug`** — remove `stale(out)`, unused `elem_out_*` accumulators

---

## Build & verify commands

```bash
# Build (CUDA release)
cmake --preset gnu-release-cuda
cmake --build --preset build-gnu-release-cuda --parallel

# Gtests
./build/gnu-release-cuda/bin/maxwell_solver_tests --device cuda \
  --gtest_filter='HesthavenOperatorBenchmark.GPU_Mult*'

# End-to-end case (watch out debug lines)
./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_RCS_PEC_Circle_G1_1m_Uniform_Hesthaven/2D_RCS_PEC_Circle_G1_1m_Uniform_Hesthaven.json \
  --device cuda 2>&1 | rg 'out  debug|TFSF debug'

# CPU baseline (same JSON)
./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_RCS_PEC_Circle_G1_1m_Uniform_Hesthaven/2D_RCS_PEC_Circle_G1_1m_Uniform_Hesthaven.json \
  --device cpu 2>&1 | rg 'out  debug|TFSF debug'
```

---

## Key code locations

| Area | Location |
|------|----------|
| GPU mult kernel | `src/evolution/HesthavenEvolution.cu` — `hesthaven_mult_gpu` |
| Mult routing | `src/evolution/HesthavenEvolution.cpp` — `Mult()`, `MultGPU()`, `MultCPU()` |
| CPU reference loop | `MultCPU()` ~L861–904 (element flux, LIFT, scatter) |
| GPU init / sync | `initGPUData()`, `hesthaven_sync_gpu_static_data()` |
| Global sync reference | `src/evolution/GlobalEvolution.cpp` ~L1169–1177 (`out.ReadWrite()`) |
| Benchmark tests | `test/maxwell/solver/HesthavenOperatorBenchmarkTest.cpp` — `GPU_Mult_*` |
| ODE state device flag | `src/evolution/Fields.h` — `all_dofs_.UseDevice(true)` |

---

## Debug output to expect today (CUDA, broken)

```
TFSF debug | abs(E)=... abs(H)=...   # grows over time
out  debug | abs(out)=0.000e+00 max(out)=0.000e+00 stale(out)=0.000e+00 ...
```

## Debug output target (CUDA, fixed)

```
out  debug | abs(out)=...              # nonzero, same order as CPU early steps
```

---

## Constraints to remember

- Do **not** edit the plan file
- Do **not** create alternate case JSONs; folder/json/mesh names must match
- Keep TFSF on **CPU** (`evaluateTFSF` on host jumps + `Write()`)
- Hesthaven scope only unless verification proves Solver export sync is still broken
