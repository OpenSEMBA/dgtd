# Session — Centered flux, upwind_alpha, and derivative split (2026-05-28)

> **Status:** Derivative split **reverted** (2026-05-27). PML ψ driver upwind blocks added; post-revert runs in [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md).

Follow-up to [12-session-pml-boundary-terminator.md](./12-session-pml-boundary-terminator.md).

## Hesthaven flux split (dgtd convention)

| Integrator | Interior coefficient | Role |
|------------|---------------------|------|
| **Directional** (volume) | 1 | Always in curl |
| **OneNormal** | **1.0** (fixed) | Centered face flux — always present |
| **ZeroNormal** | `upwind_alpha` | Upwind stabilization |
| **TwoNormal** | `upwind_alpha` | Upwind stabilization |

At **`upwind_alpha = 0`**, global cross-curl in PML uses **volume + OneNormal** only — same structure as the PML ψ driver (volume `DerivativeIntegrator` + σ-weighted OneNormal SBP).

At **`upwind_alpha = 1`**, extra Zero/Two terms act on field self-blocks; the ψ driver never included them → flux mismatch (see prior analysis).

## New test cases

| ID | Path | `upwind_alpha` |
|----|------|----------------|
| A0 | [`1D_PML_centered/`](../../testData/maxwellInputs/1D_PML_centered/) | 0.0 |
| A1 | [`1D_PML/`](../../testData/maxwellInputs/1D_PML/) | 1.0 (baseline) |
| D0 | [`1D_PML_buffer_centered/`](../../testData/maxwellInputs/1D_PML_buffer_centered/) | 0.0 |
| D1 | [`1D_PML_buffer/`](../../testData/maxwellInputs/1D_PML_buffer/) | 1.0 |

Point probes at **x = 1.99** (vacuum), **2.01** (PML, past interface), **2.99** (deep PML / near old outer SMA).

## Phase 1 — Diagnostics (no deriv split, `PML_DERIV_SPLIT` default was N/A)

`max|psi|` at Mult call **80000** (t = 20, dt = 0.001):

| Case | `upwind_alpha` | max \|psi\| @ t=20 | Notes |
|------|----------------|-------------------|--------|
| A0 | 0.0 | **6.38×10⁴** | Ringing / ψ growth (user: fields look OK mid-run, late spread) |
| A1 | 1.0 | **1.77×10¹³** | Same as prior sign audit |
| D0 | 0.0 | **9.79×10³⁰** | Buffer + centered still catastrophic |
| D1 | 1.0 | **1.25×10¹⁵** | Buffer + upwind |

**Takeaway:** Centered flux alone does **not** fix late-time ψ on buffer; it changes the failure mode vs upwind baseline.

## Phase 2 — Derivative split (reverted)

Implemented then **removed** — see [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md).

| Finding | Detail |
|---------|--------|
| False success | `max|psi| ≡ 0` with split on; vacuum **Ey** blew up immediately at the interface. |
| Root cause | Marker subtraction removed face blocks that still coupled **vacuum DOFs** at x = 2. |
| Replacement | No `PML_DERIV_SPLIT`; global curl restored on full domain. ψ driver extended with σ-weighted Zero/Two when `upwind_alpha > 0`. |

## Product guidance (updated)

1. **`upwind_alpha: 0.0`** — still best late ψ behavior in 1D matrix; does not pass probe gate at t = 20.
2. **Do not use marker-based derivative split** without interface-safe assembly.
3. Outer BC guidance unchanged from [12-session-pml-boundary-terminator.md](./12-session-pml-boundary-terminator.md).

## Recommended next steps

See [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md).
