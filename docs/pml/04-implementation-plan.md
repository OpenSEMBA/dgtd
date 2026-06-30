# Implementation plan and milestones

## Git workflow

1. Checkout **`dev`** (or user-specified base branch).
2. **Delete remote/local `vol_pml`** if it exists (user requested wipe).
3. Create fresh **`vol_pml`** from current dev.
4. All Milestone A work lands on **`vol_pml`** via focused commits (user must request commits explicitly per project rules).

**Do not** cherry-pick from old `vol_pml` implementation.

---

## Milestone A — Generic CFS-CPML on CPU serial (RK4)

**Definition of done:** User's 1D mesh runs; vacuum probes show reflected field **≤ −40 dB** vs incident (user DFT); no regression when no PML tags; `SBC_PML` removed.

### Step A.1 — Cleanup (no new physics)

| Task | Details |
|------|---------|
| Remove `SBC_PML` | `driver.cpp`: `isSGBCBoundaryType`, `buildAutoPMLSGBCLayers`, boundary parsing |
| Remove README `SBC_PML` docs | Replace with pointer to `docs/pml/` |
| Audit dead PML code | `VolumetricRegionSubMesher` may remain unused; remove only if orphaned imports/tests break |
| Clean test JSON | Remove or update legacy boundary PML cases (`1D_SGBC_PML`, `*_G1_PML`, `*_G2_PML`) per user — **volumetric only** going forward |
| Keep | `2D_RCS_Circle_Vol_PML` as structural reference (update when parser lands) |

**Verify:** Project builds; existing non-PML cases still run.

### Step A.2 — Model and JSON

| Task | Details |
|------|---------|
| Add `PMLProperties.h/cpp` | Struct + axis parsing |
| Extend `Model` | `hasPML()`, region list, tag → property lookup |
| Extend `driver.cpp` | `"type":"PML"`, `"type":"vacuum"` |
| Validation | Parse `2D_RCS_Circle_Vol_PML.json` without error |

**Verify:** Driver loads PML tags; no simulation change yet.

### Step A.3 — Geometry tables at init

| Task | Details |
|------|---------|
| PML element mask | From attributes |
| QP profile tables | κ, σ, α per (element, qp, direction) |
| Depth ρ_d | Neighbor/geometry algorithm ([03-json-schema-and-mesh.md](./03-json-schema-and-mesh.md)) |
| Interface check | Log max \|σ\|, \|α\|, \|κ−1\| at interface QPs → should be ~0 |

**Verify:** Diagnostic print on rank 0 for reference case; interface values near zero.

### Step A.4 — Extended state

| Task | Details |
|------|---------|
| `PMLAuxLayout` | Compute `n_aux`, offsets |
| `Fields` | `SetSize(6*N + n_aux)`; ψ block after E/H |
| `GlobalEvolution` ctor | `TimeDependentOperator(6*N + n_aux)` |
| `Mult` | Assert/resize `in`, `out` to extended size |
| Probes | Confirm only touch first `6*N` |

**Verify:** No-PML case: `n_aux=0`, identical to previous behavior.

### Step A.5 — Operators (physics)

| Task | Details | Status |
|------|---------|--------|
| κ in PML mass | `buildEpsMuPiecewiseVector` or equivalent | Deferred (Milestone A uses `kappa_max=1`) |
| `PMLOperator_` / `applyPMLCoupling` | Factory-built ADE CFS-CPML via `DGOperatorFactory::buildPMLOperator` | **Done** |
| Wire in `Mult()` | After `globalOperator_->Mult`, before TFSF | **Done** |
| Initialize ψ = 0 | At simulation start | **Done** |

**Verify:** Simulation runs without NaN on first RK stages; ψ activates when the pulse enters PML. Long runs may need reduced `time_step` (explicit stiffness from σ); see Milestone B ATS/SDIRK.

### Step A.5b — Factory `PMLOperator_` blocks

| Block | Factory builder | Status |
|-------|-----------------|--------|
| ψ mass \(-\alpha\) | Per-component `ψ^E_{d,c}` / `ψ^H_{d,c}` blocks: `PMLProfileCoefficient::Alpha` + marked mass + scalar `M^{-1}` | **Done** (zero when `alpha_max=0`) |
| ψ driver \(\sigma \mathcal{D}_d\) | Per active component: marked volume + face integrators, driver column from `pmlPsiEDriverCoupling` / `pmlPsiHDriverCoupling`; SBP volume/face opposite sign on same column | **Done** |
| Field correction \(\pm \psi/\kappa\) | One field row per component block (`H_c ← ψ^E_{d,c}`, `E_c ← ψ^H_{d,c}`); Maxwell `M^{-1}[E/H]` on field rows | **Done** |

### Step A.6 — User acceptance

| Task | Details | Status |
|------|---------|--------|
| User 1D case | `testData/maxwellInputs/1D_PML/` | Ready |
| Probes in vacuum | x = −2.75, 1.0 | Ready |
| DFT | User offline; target **−40 dB** | **Pending user sign-off** |
| `1D_PML` boundary tag 2 | SMA restored for acceptance (was PEC for bounce smoke test) | **Done** |

**Verify:** User/Salvador sign-off on reflection level at `time_step ≤ 0.001` for the reference case.

---

## Milestone B — Hardening (post-A)

| Item | Notes |
|------|-------|
| `ImplicitSolve` + SDIRK | Include PML blocks in `(I - dt*A)` when `alpha_max` stiff |
| `2D_RCS_Circle_Vol_PML` validation | Same −40 dB methodology |
| ATS adjustment | CFL may need PML-related limit on α |
| README complete | All JSON fields documented |

---

## Milestone C — Production extras (deferred)

| Item | Notes |
|------|-------|
| GPU `GlobalEvolution.cu` | Extended load/scatter kernels |
| MPI regression | Extended vector ghost exchange if ψ lives on same FES |
| Automated tests | `test/cases` if desired later |
| JSON schema file | `core/schema` if project adds it |

---

## Implementation order (strict)

```text
A.1 cleanup → A.2 parser → A.3 profiles → A.4 state → A.5 operators → A.6 user test
```

Do **not** implement operators before extended state (Step A.4 before A.5).

Do **not** implement GPU before CPU validation.

---

## Risk register

| Risk | Mitigation |
|------|------------|
| Wrong ADE signs | Transcribe Gedney & Zhao Maxwell equations; compare 1D limit to literature |
| Interface reflection | Enforce σ=0, κ=1, α=0 at ρ=0; check vacuum neighbor detection |
| RK4 instability with α | Reduce dt; Milestone B SDIRK |
| SMA outer termination + α=0 late-time ψ drift | Localized HF noise at outer PML edge (x=3 in 1D); try `alpha_max > 0`, smaller dt, or SDIRK (Milestone B) |
| Hard-coded `6*ndofs` | Grep and fix all locations listed in [02-codebase-architecture.md](./02-codebase-architecture.md) |
| Corner double-counting | Independent ρ_X, ρ_Y; full Maxwell ADE corners |
| Confusion with bulk σ | Separate code paths; static_assert PML tags have zero bulk σ |

---

## Agent session strategy

Each Agent session should:

1. Read `docs/pml/README.md` and relevant sub-doc.
2. State which step (A.1–A.6) is in scope.
3. Stop after step verification criteria met.
4. Update this file's checkboxes in **Open items** sections of other docs when resolved.

Suggested first Agent prompt after docs:

> Implement Step A.1 + A.2 on fresh `vol_pml` from dev per `docs/pml/04-implementation-plan.md`.

---

## Open checklist (update as work proceeds)

- [x] A.1 Cleanup complete
- [x] A.2 Parser complete
- [x] A.3 Profile tables complete
- [x] A.4 Extended state complete
- [x] A.5 PMLOperator / ADE complete
- [ ] A.6 User −40 dB acceptance (**partial:** bulk absorption good to t≈12; **t=20 unstable** for α=0 and α=1 under RK4; probe/DFT workflow updated — see [09-session-record-2026-05-27.md](./09-session-record-2026-05-27.md), [10-session-2026-05-28-stability-dft.md](./10-session-2026-05-28-stability-dft.md))
- [ ] Gedney & Zhao equation numbers documented in `01-physics-and-formulation.md`
- [ ] σ_max formula validated
- [ ] README updated
