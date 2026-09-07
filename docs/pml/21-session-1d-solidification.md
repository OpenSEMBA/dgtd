# Session record — 2026-09-04 (1D solidification)

Follow-up to [`20-pipeline-reexploration.md`](./20-pipeline-reexploration.md). Commits on `dev`:

| Commit | Content |
|--------|---------|
| `bf37b094` | ADE pole \(\alpha+\sigma/\kappa\), field-correction signs, DFT tooling, L=1/2/3 + SMA cases |
| `0b4ef20f` | Disable σ-weighted Zero/Two in ψ ADE driver |

## Locked 1D picture (do not re-litigate without new evidence)

### Stability (late-time)

Fixed by Gedney ADE transcription:

1. ψ mass = **Decay** \(\alpha + \sigma/\kappa\) (not \(\alpha\) alone).
2. Field corrections: \(\dot E \mathrel{-}= \psi^H/\kappa\), \(\dot H \mathrel{+}= \psi^E/\kappa\).

`1D_PML_buffer` / `1D_PML` remain bounded through at least \(t=60\) with post-pulse \(\max|\psi|\) at noise.

### Acceptance (−40 dB DFT)

Tooling: [`scripts/pml_dft_reflection.py`](../../scripts/pml_dft_reflection.py).

Primary case: [`1D_PML_DFT`](../../testData/maxwellInputs/1D_PML_DFT/) (L=1, probe at \(x=0\), inc [3.0, 7.2], ref [8.5, 11.5]).

| Configuration | \(R_{\mathrm{dB}}(f_{\mathrm{peak}})\) | Gate |
|---------------|----------------------------------------|------|
| ADE fixed + **ψ Zero/Two on** (global upwind) | ≈ −29 | FAIL |
| ADE fixed + **ψ Zero/Two off** (global upwind kept) | ≪ −40 (~−300, noise) | **PASS** |
| Same L, `upwind_alpha=0` | ≈ −87 | PASS |
| Same L, order 3, ψ-upwind on | ≈ −40 | near miss |
| Same L, dx=0.05 (20 els), ψ-upwind on | ≈ −135 | PASS |
| L=2 / L=3 (any m, R in sweep), either ψ-upwind mode | ≪ −40 | PASS |
| SMA-only twin [`1D_SMA_DFT`](../../testData/maxwellInputs/1D_SMA_DFT/) | noise floor | PASS (1D normal incidence exact) |

**Production code choice:** `pml_upwind = false` in `collectPMLOperatorBlocks` (`DGOperatorFactory.h`). Global `upwind_alpha` is unchanged.

## What “ψ upwind off” actually means

**Not** “PML has no upwind.” Maxwell still uses Hesthaven upwind in `globalOperator_` wherever `upwind_alpha > 0`.

| Operator | Volume \(D\) | OneNormal (SBP) | Zero/Two upwind |
|----------|--------------|-----------------|-----------------|
| `globalOperator_` (Ė, Ḣ) | yes | yes | yes if `upwind_alpha>0` |
| `PMLOperator_` ψ driver | yes (σ/κ) | yes (σ/κ, oppose vol) | **no** (disabled) |
| `PMLOperator_` field corr | — | — | ±ψ/κ only |

Continuous ADE wants ψ driven by the **same** discrete \(D\) as the curl. Separately assembling σ-weighted Zero/Two with a **PML volume marker** does **not** reproduce that upwind part cleanly on vacuum–PML faces; it corrupted interface matching on coarse/thin layers.

Omitting ψ-upwind leaves a **mild inconsistency** (ψ tracks centered-SBP \(D\); fields use upwind curl) that is empirically excellent in 1D. A future “proper” fix is to share the global discrete \(D\) blocks into the ψ driver (not to re-enable the old marker-based Zero/Two path).

## Thickness vs resolution (important for mesh design)

Early L=1 FAIL vs L=2 PASS looked like “must pre-cook a thick PML.” Ablations show:

- Controlling factor was largely **elements / profile resolution through the layer** and **broken ψ-face upwind**, not σ budget.
- On L=1 with ψ-upwind still on, **smaller `target_reflection` (larger σ) made DFT worse**.
- After disabling ψ-upwind, **L=1 passes**; thickening is optional margin, not a hard requirement for −40 dB on this pulse.

Volumetric \(L\) remains a real continuous parameter (Taflove-style), but the cliff was mostly discrete.

## SMA comparison (1D only)

[`1D_SMA_DFT`](../../testData/maxwellInputs/1D_SMA_DFT/) with two SMAs and no PML is essentially perfect for **1D normal-incidence TE** (Silver–Müller exact). That baseline is **not** a fair contest for judging PML value; use it only as a numerical floor. PML must prove itself in **2D / oblique** next.

## Case map (1D)

| Path | Role |
|------|------|
| `1D_PML_DFT` | Primary −40 dB DFT (L=1) |
| `1D_PML_L2` / `1D_PML_L3` | Thick-L confirmation |
| `1D_SMA_DFT` | SMA-only floor (1D) |
| `1D_PML_buffer` | Stability / Mult-diag workhorse |
| `1D_PML`, `*_centered`, `*_PEC`, `*_PML_NONE` | Earlier matrix (still useful) |

Optional local experiment mesh `1D_PML_L1fine` (dx=0.05) may exist untracked — same L, finer PML spacing.

## Open for 2D (next session — banter, not locked)

1. Pick a deciding 2D geometry (oblique / grazing where SMA is weak).
2. Decide probe + DFT protocol in 2D (or power flux / RCS).
3. Whether shared-\(D\) ψ driver is needed before or after first 2D −40 dB attempt.
4. Corner tags (`active_axes` X+Y) already in `2D_RCS_Circle_Vol_PML` JSON shape.

## Do not retry without new theory

- Re-enable marker-based ψ Zero/Two “to restore upwind consistency”
- Blind σ / `target_reflection` increases on under-resolved L=1 as a reflection fix
- Marker `PML_DERIV_SPLIT` / regional curl omit (already reverted historically)
- Using 1D SMA-vs-PML as the quality bar for volumetric PML
