# Verification and acceptance

## Success criterion (locked)

**−40 dB** reflected power relative to incident, evaluated in **frequency domain** from **vacuum probe** time series (user performs DFT offline).

This is slightly looser than JSON `target_reflection: 1e-6` (−60 dB design target) — acceptance uses **−40 dB** for Milestone A pass.

---

## Probe strategy

### Placement

- Probes **only in vacuum region** (non-PML element attributes).

**1D `1D_PML` reference layout (corrected 2026-05-28):**

| Probe | x | Role |
|-------|---|------|
| PointProbe0 | **−2.75** | Intended **reflection** (upstream of TFSF at x = −2.5) — see note below |
| PointProbe1 | **1.0** | **Incidence**, then **reflection** (separate by time window) |

- **|E_inc|:** Ey on PointProbe1, early window **t ∈ [3.5, 6.5] s** (normalized time; use step index × `time_step` if parsing `.dat` time column).
- **|E_ref|:** Ey on PointProbe0 when nonzero, **or** late window on PointProbe1 (e.g. t ≳ 8 s). Incidence at x = 1.0 ends ~t ≈ 6.5 s before the reflected lobe grows.

**Note:** Exports to date show **Ey ≡ 0** at x = −2.75 (upstream of TFSF); use PointProbe1 late-time window for |E_ref| until that is resolved. Details: [10-session-2026-05-28-stability-dft.md](./10-session-2026-05-28-stability-dft.md).

User responsibility: probe coordinates in JSON for their 1D Gmsh case.

### Fields recorded

- **E and/or H components** as supported by existing `PointProbe` / exporter.
- **Not** auxiliary ψ.

### Post-processing (script)

Use [`scripts/pml_dft_reflection.py`](../../scripts/pml_dft_reflection.py) on vacuum PointProbe
exports. Dedicated case: [`testData/maxwellInputs/1D_PML_DFT/`](../../testData/maxwellInputs/1D_PML_DFT/).

1. Run simulation to `final_time` sufficient for the reflected lobe to clear the probe.
2. Export probe time series (existing probe/export pipeline).
3. DFT incident and reflected **time windows** separately (Hann taper default).
4. Compute \(R(f) = |E_{\mathrm{ref}}(f)| / |E_{\mathrm{inc}}(f)|\).
5. **20 log10(|R|) ≤ −40 dB** at the incident spectral peak → pass.

```sh
mpirun -np 1 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i testData/maxwellInputs/1D_PML_DFT/1D_PML_DFT.json
python3 scripts/pml_dft_reflection.py Exports/single-core/1D_PML_DFT --probe 0
```

**2026-09-04 (solidified):** After ADE pole/sign fixes, L=1 `1D_PML_DFT` **FAIL**ed (~−29 dB)
with σ-weighted Zero/Two in the ψ driver, and **PASS**es (≪ −40 dB) with those ψ-upwind
blocks disabled. Global `upwind_alpha` is unchanged — Maxwell remains upwind; only the ADE
driver drops marker-based Zero/Two. See [`21-session-1d-solidification.md`](./21-session-1d-solidification.md).

Related cases: [`1D_PML_L2`](../../testData/maxwellInputs/1D_PML_L2/), [`1D_PML_L3`](../../testData/maxwellInputs/1D_PML_L3/),
[`1D_SMA_DFT`](../../testData/maxwellInputs/1D_SMA_DFT/) (1D SMA floor; not a fair PML contest).

---

## Test cases

### Primary acceptance (user-provided)

| Case | Path | Purpose |
|------|------|---------|
| 1D PML slab | `testData/maxwellInputs/1D_PML/` | First −40 dB check |
| Mesh | Gmsh 1D: vacuum + right-side PML tag | Minimal debugging |

Suggested setup:

- Pulse source: TFSF planewave (`1D_PML.json`: Ey, +x, tag 3 at x = −2.5).
- PML on **+X** side only: vol tag 3, `"active_axes": ["X"]`.
- Outer **left** (x = −3): SMA; outer **right** (x = 3): SMA.
- Probes in vacuum only: x = −2.75 (reflection target), x = 1.0 (incidence + reflection by time window).

### Structural reference (not necessarily passing until implemented)

| Case | Path |
|------|------|
| 2D vol PML RCS | `testData/maxwellInputs/2D_RCS_Circle_Vol_PML/` |

Use for **JSON shape**, multi-block `active_axes`, corner tags — Milestone B validation.

### Regression (no PML)

Run existing cases without PML tags after each milestone:

- e.g. `testData/maxwellInputs/1D_PEC/`
- Confirm **bit-identical or numerical match** within tolerance when `n_aux = 0`.

---

## Unit-level checks (developer, no automated CI required)

### Interface profile check (init)

For each PML region, sample QPs nearest vacuum neighbor:

```text
|κ - 1| < tol
|σ| < tol
|α| < tol
```

tol suggestion: 1e-12 normalized.

### Zero PML overhead

Model with **no** `"type":"PML"`:

```text
n_aux == 0
TimeDependentOperator::Height() == 6 * ndofs
Mult() behavior unchanged (regression)
```

### Energy sanity (qualitative)

During 1D run:

- No exponential blow-up in first few periods.
- Fields in PML decay in amplitude with depth (visual in Paraview on E/H).

### NaN / Inf guard

After N time steps, `fields.allDOFs().Norml2()` finite.

---

## Comparison baselines (optional)

| Baseline | Use |
|----------|-----|
| Larger PML thickness in Gmsh | Reflection should not worsen when thickness adequate |
| `alpha_max = 0` vs `> 0` | Long-time: CFS should reduce late-time drift |
| Remove PML tag (SMA only) | Should show stronger boundary artifacts — confirms PML helps |

**Not** comparing to old `SBC_PML` (removed).

---

## Failure triage guide

| Symptom | Likely cause |
|---------|--------------|
| ~0 dB reflection | Interface σ≠0 or wrong neighbor; PML not applied |
| Slow linear growth | Wrong sign in ADE; RK4 dt too large |
| Late-time blow-up | α=0 classical PML; need CFS or SDIRK |
| Direction wrong | `active_axes` mismatch with mesh orientation |
| Corner spikes | Corner tag missing one axis profile |
| Good 1D, bad 2D | Corner ADE terms; independent ρ per axis |

---

## Documentation of results

When Milestone A passes, record in PR or `docs/pml/05-verification.md` appendix:

- Case name, mesh file, dt, order, α_max, κ_max, grading_order
- Probe locations
- Measured dB reflection at frequency f
- Date and commit hash

---

## Explicit non-requirements

- No automated gtest in `test/cases/` for Milestone A
- No Paraview images in repo required
- No comparison to analytical Green's function unless user adds

---

## Acceptance sign-off

| Role | Criterion |
|------|-----------|
| User | −40 dB on 1D case with their DFT |
| Salvador (optional) | Formulation matches Gedney & Zhao ADE CFS-PML intent |
| Agent/dev | Regression on no-PML cases; no SBC_PML remains |
