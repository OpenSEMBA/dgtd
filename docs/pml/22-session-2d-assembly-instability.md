# Session record — 2026-09-04 (2D assembly instability)

Follow-up to [`21-session-1d-solidification.md`](./21-session-1d-solidification.md).

## Cases

| Case | Role |
|------|------|
| `2D_PML_X_slab` | X-only slab from `2D_InteriorBoundary_SGBC_Test` mesh (TM `Ez`) |
| `2D_PML_X_slab_centered` / `_TE` / `_short` | Flux / polarization / window variants |
| `2D_RCS_Circle_Vol_PML_*` | Multi-axis circle (still secondary) |

## Ablation map (`2D_PML_X_slab_short`, before nodal decay)

| Parts | Fields | ψ |
|-------|--------|---|
| `PML_SKIP_PML` | O(1) stable | (no diag) |
| driver only | O(1) | peak ~85, bounded |
| mass + driver (no corr) | O(1) | **explodes** |
| mass + corr (no driver) | O(1) | 0 |
| full | explodes | explodes |

**Conclusion:** open-loop `M⁻¹ Mass(σ)` decay is unstable in 2D; closed-loop corrections then wreck the fields.

## Spectral evidence

Mass-only ψ block eigenvalues (`PML_SKIP_DRIVER=1 PML_SKIP_CORR=1`, export CSR):

| Mesh | eig min | eig max | # positive |
|------|---------|---------|------------|
| 1D `1D_PML_buffer` | −35.3 | **+4e−6** | ~8 (noise) |
| 2D X-slab | −39.3 | **+18.5** | **27** |

So `−M⁻¹ Mass(σ)` is **not negative-semidefinite** on 2D triangles. That is the ψ-only blow-up.

## Face-marker API bug (fixed, not the stabilizer)

`AddInteriorFaceIntegrator(integ, marker)` in mfem-geg treats `marker` as a **boundary-attribute ignore** list. PML passed a **volume** `pml_marker` (size mismatch on 2D). Face integrators now use `buildInteriorIgnoreMarker()` like global; σ=0 outside PML still localizes.

## Fixes landed

1. Face marker semantics (hygiene).
2. **dim ≥ 2:** nodal Gauss–Lobatto diagonal for ADE decay `(α+σ/κ)` and for field-correction `1/κ`.
3. **dim == 1:** keep validated `M⁻¹ Mass(σ)` / Maxwell `M⁻¹ Mass(1/κ)`.
4. Richer Mult diag (`Ez`, per-component `psiE/H`).
5. Ablation envs: `PML_SKIP_MASS`, `PML_SKIP_DRIVER`, `PML_SKIP_CORR`, `PML_SKIP_FACE_DRIVER`.

After (2): **mass+driver is stable** on 2D X-slab (peak ψ ~15, Ez O(1)).

## Still open

**Full operator (driver + nodal decay + corrections) still grows** on 2D X-slab (secular / explosive depending on signs). `PML_SIGN_TEST=1` (flipped corrections) slows growth but is not a fix and hurts 1D absorption if applied there.

Likely next: make the **driver** consistent with nodal decay (same discrete `D`), or shared global `D` blocks into ψ — not more σ tuning.

## Do not

- Re-enable marker-based ψ Zero/Two on 2D without new evidence
- Use under-integrated `InverseIntegrator` IR (singular mass)
- Treat 1D SMA vs PML as the 2D quality bar
