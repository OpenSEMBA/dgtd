# Session — Centered case: frame ↔ probe correlation (2026-05-27)

Cross-check of user ParaView frames (`cent.0000`–`cent.0199`, 30 shots) against point probes and `PML_OPERATOR_AUDIT` for [`1D_PML_centered`](../../testData/maxwellInputs/1D_PML_centered/1D_PML_centered.json).

## Time mapping (important)

| Quantity | Value |
|----------|--------|
| JSON `final_time` | **20** (normalized code units, c = 1) |
| JSON `time_step` | **0.001** → **20 000** steps |
| Displayed in stats / probes | **66.713 ns** (= `20 / c_SI × 10⁹`) |
| ParaView cycles | **0 … 199** (`exporter.saves: 200`) |
| **t_code(cycle k)** | `20 × k / 199` |
| **t_ns(cycle k)** | `66.713 × k / 199` |

Probe `.dat` time column is **`t_code / c_SI`** (meters), not seconds.

| User frame | Cycle | t_code | t_ns |
|------------|-------|--------|------|
| Vacuum pulse “just right” | **48**, **61** | 4.82, 6.13 | **16.1**, **20.5** |
| Peak near **x = 2.0** | **68 – 75** | 6.83 – 7.54 | **22.8 – 25.1** |
| Late blow-up | **199** | 20.0 | **66.7** |

## Domain (mesh x → plot x ≈ mesh x + 3)

| Tag | Region | Mesh x |
|-----|--------|--------|
| 1 | Left buffer | −3 … −2.5 |
| 2 | Vacuum | −2.5 … **2** |
| 3 | PML | **2** … 3 |

Probes: **1.99** (vacuum), **2.01** (PML), **2.99** (deep PML). Source: Gaussian planewave in tag **3**, +x propagation (energy enters vacuum at x = 2).

## Probe amplitudes at user frames

`|Ey|` and `|Hz|` from `Exports/single-core/1D_PML_centered/PointProbes/` (TE: dominant component at probes).

| Cycle | t_ns | Probe0 (1.99) | Probe1 (2.01) | Probe2 (2.99) |
|-------|------|---------------|---------------|---------------|
| 48 | 16.1 | ~0 | ~0 | ~0 |
| 61 | 20.5 | **0.014** | **0.012** | ~10⁻⁵ |
| 68 | 22.8 | **0.50** | **0.47** | **0.016** |
| 75 | 25.1 | **0.85** | **0.87** | **0.84** |
| 102 | 34.2 | ~0* | ~0* | ~0.02 |
| 199 | 66.7 | **~513** | **~285** | **~1793** |

\*Cycle 102 is near a zero crossing in the short probe record, not necessarily quiescent globally.

### Interpretation vs your visual picks

1. **Frames 48 / 61 (vacuum-travelling pulse)** — Line plots can show a smooth lobe in tag 2 while **1.99 stays small** early: the packet is still **right of the probe** (source at x > 2, wave propagates −x). By **61**, probes are O(10⁻²–10⁻²) — consistent with a **small but arriving** signal, not yet the interface peak.

2. **Frames 68 – 75 (maximum near x = 2)** — Probes bracket the interface: **|Ey| ≈ 0.47–0.87** at **2.01** and **1.99**, peak in the **68–75** window. This matches “maximum at 2.0” in the **O(1) physical** sense, **before** late instability.

3. **Instability is later** — First probe sample with **|Ey| or |Hz| > 3** (gate threshold): **cycle ≈ 167** (t_ns ≈ **56**), not at 68–75. So the interface window is **not** the blow-up instant; it is the last **clean** interface-dominated epoch.

4. **`check_pml_probe_stability.py`** — **FAIL** on this export (max |Ey|, |Hz| at probes **~7×10²–2×10³**), dominated by late-time growth, not the 68–75 window.

## Assembly audit (post-revert, centered α = 0)

`PML_OPERATOR_AUDIT=1` on init:

| Metric | Value |
|--------|--------|
| Global nnz | 3360 |
| `||Hz←Ey||_F` (full) | 776.46 |
| `||Hz←Ey||_F` on PML DOFs | **314.55** |
| `||ψ←Ey||_F` | 3677.07 |
| `driver_upwind_nnz` | **0** (α = 0) |
| PMLOperator nnz | 936 |

Duplicate-curl architecture (full global curl + σ ψ driver) is unchanged after reverting regional assembly.

## Conclusions

| # | Finding |
|---|---------|
| 1 | User frame labels are **consistent with probe data** for the **interface peak (68–75)**. |
| 2 | Vacuum frames **48/61** are **earlier** than interface peak; low amplitudes at x = 1.99 do **not** contradict a visible spatial pulse farther right in vacuum. |
| 3 | **Stable O(1) physics** appears limited to roughly **t_ns ≲ 56** (cycle ≲ 167); debugging should prioritize **ψ signs / coupling** and additive PML formulation, not ParaView export. |
| 4 | Do not read JSON `final_time: 20` as 20 seconds; displayed run length is **~67 ns**. |

## Suggested next checks (not run here)

- Export **ψ** or `max|psi|` vs cycle alongside these frames (cycle **167** onset).
- Compare same table for **`1D_PML`** (α = 1) at cycles **48 / 61 / 68**.
- Sign sweep (session 11) focused on interface peak window **68–75** where fields are still trustworthy.
