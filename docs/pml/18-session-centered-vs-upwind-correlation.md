# Session — Centered vs upwind: frame window + ψ diagnostics (2026-05-27)

Follow-up to [17-session-centered-frame-probe-correlation.md](./17-session-centered-frame-probe-correlation.md). Fresh runs with `[PML Mult diag]` captured in `pml_diag.log` per case.

Tool: [`analyze_pml_frame_probes.py`](../../testData/maxwellInputs/1D_PML/analyze_pml_frame_probes.py)

## User frame anchors (confirmed)

| Role | Cycles | t_ns | Mult call (≈) |
|------|--------|------|----------------|
| Vacuum-travelling pulse | **48**, **61** | 16.1, 20.5 | 19300, 24500 |
| Peak near **x = 2** | **68 – 75** | 22.8 – 25.1 | 27300 – 30200 |

Probe time column = `t_code / c_SI`; `t_ns = 66.713 × cycle / 199`.

---

## Probe amplitudes at anchor frames

| Cycle | Centered P0 (1.99) | Centered P1 (2.01) | Upwind P0 | Upwind P1 |
|-------|-------------------|-------------------|-----------|-----------|
| 48 | ~0 | ~0 | ~0 | ~0 |
| 61 | 0.014 | 0.012 | 0.014 | 0.012 |
| **68** | **0.50** | **0.47** | **0.50** | **0.47** |
| **75** | **0.85** | **0.87** | **0.85** | **0.87** |
| 102 | ~0* | ~0* | 5.6 | 8.1 |
| 199 | ~513 | ~285 | ~10²⁴ | ~10²⁴ |

\*Zero crossing in probe record, not necessarily global quiescence.

**Takeaway:** Through **cycle 75**, centered and upwind are **indistinguishable** at the three probes. Divergence appears **after** the interface-peak window.

---

## ψ diagnostics at anchor frames (`[PML Mult diag]`)

| Cycle | Mult log | Centered `max|psi|` | Centered `outer_Ey` | Upwind `max|psi|` | Upwind `max|dpsi/dt|` |
|-------|----------|---------------------|---------------------|-------------------|----------------------|
| 48 | 19500 | 1.7×10⁻⁶ | 5×10⁻⁸ | 1.1×10⁻⁹ | — |
| 61 | 24500 | 3.5×10⁻⁴ | 1×10⁻⁵ | 3.5×10⁻⁴ | — |
| **68** | 27500 | **0.15** | **0.02** | **0.15** | — |
| **75** | 30000 | **2.41** | **0.72** | **2.41** | **0.50** |
| 102 | 41000 | 0.10 | 0.016 | **6.2×10³** | — |
| 144 | 58000 | 4.0 | 0.13 | **1.7×10¹⁴** | — |
| 167 | 67000 | **189** | 5.1 | **3.7×10¹⁹** | **3.0×10⁴** |
| 199 | 80000 | **6.4×10⁴** | 2530 | **2.3×10²⁷** | **5.7×10⁹** |

At **68–75**, ψ and outer-PML Ey match between cases → **interface peak is not an upwind artefact**.

---

## Instability onset (relative timing)

| Metric | Centered (α=0) | Upwind (α=1) |
|--------|----------------|--------------|
| First probe \|E/H\| > 3 | cycle **~167** (56 ns) | cycle **~101** (34 ns) |
| `max|psi|` lull | call **48000** → **0.04** | — (monotone growth after ~35000) |
| `max|psi|` blow-up | call **67000** → 189; **80000** → 6.4×10⁴ | call **36000** → 7.4; **48000** → 8.2×10⁷ |
| End `max|psi|` @ 80000 | 6.4×10⁴ | 2.3×10²⁷ |

Upwind **global Zero/Two** (6192 vs 3360 nnz) does **not** change the **68–75** snapshot, but accelerates runaway **~1.6× earlier** in probe time and **much harder** after call ~35000.

---

## Assembly audit (unchanged additive model)

Both cases after regional-split revert:

| | Centered | Upwind |
|---|----------|--------|
| Global nnz | 3360 | 6192 |
| `||Hz←Ey||_F` PML DOFs | 314.55 | 314.55 |
| `||ψ←Ey||_F` | 3677 | 3677 |
| `driver_upwind_nnz` | 0 | 432 |

---

## Recommended debugging focus

1. **Sign / coupling sweeps** (`PML_SIGN_TEST`) bracketing **cycles 68–75** only (fields still O(1)).
2. Treat **upwind** instability as **post-interface** stiffening: compare `max|dpsi/dt|` and `max|out H|` from call **35000** onward (centered still O(1) there; upwind already `dpsi` ~ 48–174).
3. Do **not** remove global flux; duplicate-curl remains a formulation/sign issue per audit, not fixed by regional omission.

---

## Artifacts

- `testData/maxwellInputs/1D_PML_centered/pml_diag.log`
- `testData/maxwellInputs/1D_PML/pml_diag.log`
