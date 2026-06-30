# Session record — 2026-05-28 (stability t=20, probes, alpha)

Follow-up to [09-session-record-2026-05-27.md](./09-session-record-2026-05-27.md). Case: [`1D_PML/`](../../testData/maxwellInputs/1D_PML/).

## Config changes

[`1D_PML.json`](../../testData/maxwellInputs/1D_PML/1D_PML.json):

- `final_time`: **20.0** (was 12.0)
- `probes.exporter.saves`: **200** (was 120)
- `alpha_max`: **0.0** (restored after sweep; user also tested α → 1.0 manually)

## Probe roles (corrected)

| Probe | x | Role |
|-------|---|------|
| PointProbe0 | −2.75 | Intended **reflection** (upstream of TFSF at −2.5) |
| PointProbe1 | 1.0 | **Incidence** then **reflection** |

### Incidence-only window at x = 1.0

Use simulation time **index × `time_step`** (probe file time column is divided by `c_SI`, not normalized time).

| Window | Purpose |
|--------|---------|
| **t ∈ [3.5, 6.5] s** | Dominant **incident** Ey peak (~1.0 at t ≈ 6.33) |
| **t > 6.5 s** | Reflected / late content begins (|Ey| ≈ 0.91 at t ≈ 6.50, then decays) |

### PointProbe0 at x = −2.75

After fresh runs to t = 20, **all components stay zero** (Ey and Hz). Likely cause: position is **upstream of the TFSF surface** (x = −2.5); the TFSF formulation may not populate the total field there for scattered/reflected waves.

**Practical DFT workaround until probe/export is fixed:** use **PointProbe1 only** — |E_inc| from early window, |E_ref| from a late window (e.g. t ∈ [8, 12] or after the main pulse).

Crude peak ratio from exported run (probe1 only): |E_ref|/|E_inc| ≈ 1.7×10⁻⁴ in t ∈ [8, 12] → **≈ −75 dB** (peak-based, not full DFT).

## Late-time stability (`final_time = 20`, `dt = 0.001`)

### `alpha_max = 0.0`

| Simulation time (approx.) | max \|ψ\| | Notes |
|---------------------------|-----------|--------|
| t ≈ 7.5 (call 30000) | **2.4** | Matches healthy t = 12 behavior |
| t ≈ 11.25 (call 45000) | 0.17 | Still bounded |
| t ≈ 13.75 (call 55000) | **1417** | Onset of growth |
| t = 20 (call 80000) | **1.77×10¹³** | **Blow-up** |

**Conclusion:** α = 0 is **not stable to t = 20**; localized outer-edge mode becomes exponential after ~t ≈ 12–14.

### `alpha_max = 1.0`

| Simulation time (approx.) | max \|ψ\| |
|---------------------------|-----------|
| t ≈ 7.5 | 2.28 |
| t ≈ 11.25 | **0.04** (better than α = 0) |
| t ≈ 13.75 | 238 |
| t = 20 | **5.34×10¹¹** |

**Conclusion:** CFS damping **delays and reduces** mid-time ψ noise (consistent with user visual tests) but **does not prevent** late-time blow-up to t = 20 under explicit RK4.

## User observations (aligned)

- Increasing α toward 1.0 **delays / slightly reduces** SMA-edge noise vs α = 0.
- Noise **remains** compared with pure SMA absorption.
- PML acceptance (−40 dB reflection) is separate from “SMA-like” nulling at the outer boundary.

## Regression

`1D_PEC` completes with no PML extended state.

## Open next steps

1. **DFT sign-off** — user DFT with corrected windows; fix or relocate reflection probe if −2.75 must be used.
2. **Stability** — Milestone B: smaller `dt`, SDIRK with PML in `ImplicitSolve`, or stop acceptance runs at `final_time ≤ 12` until stiff integration lands.
3. **Investigate** why Ey ≡ 0 at x = −2.75 (TFSF upstream vs probe bug).
4. **Do not** treat t = 20 α = 0/1 runs as physically meaningful past blow-up onset.
