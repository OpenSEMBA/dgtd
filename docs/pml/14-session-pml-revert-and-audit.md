# Session — Revert derivative split, PML upwind driver, probe gate (2026-05-27)

Follow-up to [13-session-centered-flux-and-deriv-split.md](./13-session-centered-flux-and-deriv-split.md).

## Revert: `PML_DERIV_SPLIT` removed

Marker subtraction of stretch-direction **volume + OneNormal** from `globalOperator_` was **retracted**.

| Issue | Detail |
|-------|--------|
| Vacuum–PML interface | PML-marked face forms still couple **vacuum DOFs** at the interface (e.g. x = 2). Subtracting those blocks from the global operator **corrupts** the semidiscrete curl at the interface. |
| Symptom | `1D_PML_centered` with split on: **immediate** vacuum instability (~t ≈ 2.5); `max|psi| ≡ 0` was a false positive. |
| Code | Deleted [`PMLDerivSplit.h`](../../src/components/PMLDerivSplit.h) / [`.cpp`](../../src/components/PMLDerivSplit.cpp). `collectGlobalDirectionalOperators` and `collectGlobalOneNormalOperators` again assemble on the **full domain**. |

**Current model:** `globalOperator_` carries the full Maxwell curl in PML and vacuum; `PMLOperator_` adds σ·D_stretch + ψ mass + ψ/κ corrections (duplicate-curl concern remains — see [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md)).

## Phase 1b — σ-weighted Zero/Two in ψ driver

When `upwind_alpha > 0` (`pd_.opts.alpha`), `PMLOperator_` now adds the same upwind face flux structure as the global operator, restricted to **PML volume marker** (interior faces only):

| Block | Integrator | Global sign | PML driver |
|-------|------------|-------------|------------|
| Zero | `MaxwellDGCoefficientZeroNormalJumpIntegrator(σ × α)` | −1 on `f[c]←f[c]` | −1 on `in_field[in_c]` |
| Two | `MaxwellDGCoefficientTwoNormalJumpIntegrator(σ × α)` | +1 on `f[c]←f[c']` when first dir = c | +1 when `pair.first == in_c` |

Assembly log: `driver_upwind_nnz`, `upwind_alpha` on `[PML] PMLOperator_ assembled`.

Files: [`BilinearIntegrators.h`](../../src/mfemExtension/BilinearIntegrators.h), [`DGOperatorFactory.h`](../../src/components/DGOperatorFactory.h) (`buildPMLDomainZeroNormalSubOperator`, `buildPMLDomainTwoNormalSubOperator`, `collectPMLUpwindDriverBlocks`).

## Uniform probes and stability gate

[`1D_PML.json`](../../testData/maxwellInputs/1D_PML/1D_PML.json) and [`1D_PML_buffer.json`](../../testData/maxwellInputs/1D_PML_buffer/1D_PML_buffer.json) now use the same point probes as centered cases:

- **x = 1.99** (vacuum, upstream of interface)
- **x = 2.01** (PML, just inside interface)
- **x = 2.99** (deep PML)

Gate script: [`check_pml_probe_stability.py`](../../testData/maxwellInputs/1D_PML/check_pml_probe_stability.py) — **PASS** if max |Ey| and max |Hz| ≤ **3** on each probe over the full run.

```bash
python3 testData/maxwellInputs/1D_PML/check_pml_probe_stability.py \
  testData/maxwellInputs/1D_PML/Exports/single-core/1D_PML
```

## Run matrix (t = 20, dt = 0.001, RK4 → Mult call 80000)

Built: `build/gnu-release-mpi/bin/opensemba_dgtd`.

### `[PML] PMLOperator_` assembly

| Case | `upwind_alpha` | nnz | `driver_upwind_nnz` |
|------|----------------|-----|---------------------|
| `1D_PML` | 1.0 | 828 | 432 |
| `1D_PML_centered` | 0.0 | 936 | 0 |
| `1D_PML_buffer` | 1.0 | 846 | 480 |
| `1D_PML_buffer_centered` | 0.0 | 972 | 0 |

### `[PML Mult diag]` (selected calls)

| Case | @ 30000 | @ 48000 | @ 80000 (t=20) |
|------|---------|---------|----------------|
| `1D_PML` | 2.41 | 8.18×10⁷ | **2.26×10²⁷** |
| `1D_PML_centered` | 2.41 | **0.0405** | **6.38×10⁴** |
| `1D_PML_buffer` | 5.11 | 1.86×10¹⁴ | **1.13×10³⁹** |
| `1D_PML_buffer_centered` | 26.9 | 9.38×10¹¹ | **9.79×10³⁰** |

Centered flux shows a **mid-run lull** at call 48000 (`max|psi| ≈ 0.04`) before late ψ growth; upwind cases blow up earlier in ψ.

### Probe gate (max |Ey|, |Hz| ≤ 3)

| Case | Probe0 (1.99) | Probe1 (2.01) | Probe2 (2.99) |
|------|---------------|---------------|---------------|
| `1D_PML` | FAIL ~10²⁴ | FAIL ~10²⁴ | FAIL ~10²⁵ |
| `1D_PML_centered` | FAIL ~10³ | FAIL ~10³ | FAIL ~10³ |
| `1D_PML_buffer` | FAIL ~10³⁷ | FAIL ~10³⁷ | FAIL ~10³⁵ |
| `1D_PML_buffer_centered` | FAIL ~10²⁸ | FAIL ~10²⁸ | FAIL ~10²⁸ |

**No case passes** the t = 20 gate. Centered upwind case is orders of magnitude smaller than upwind on vacuum probe 1.99 but still far above threshold.

## Conclusions

1. **Derivative split via marker subtract is wrong** — do not reintroduce without face-local splitting that preserves vacuum–PML interface rows.
2. **PML ψ driver now includes σ-weighted Zero/Two when `upwind_alpha = 1`** — removes the documented global vs ψ upwind mismatch; **does not** stabilize t = 20 in this matrix.
3. **`upwind_alpha: 0.0`** remains the better flux choice for late ψ (centered), but **late ringing / ψ growth** persists; buffer mesh does not fix it.
4. **Duplicate curl** (`globalOperator_` + `PMLOperator_` stretch terms) and **sign/κ** audits in [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md) remain the leading hypotheses for ψ blow-up.

## Recommended next steps

1. Early-time gate: require probe 1.99 PASS for t < 7 before pulse arrival (vacuum sanity).
2. Audit σ·D driver magnitude vs global curl in PML elements (ParaView / block export).
3. Optional: regional disable of global directional+OneNormal **inside PML volume only** without subtracting interface faces (new assembly path, not marker subtract of assembled global blocks).
