# Session — Sign matrix at interface window (t = 8.05, 2026-05-27)

Short rerun of [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md) on [`1D_PML_centered`](../../testData/maxwellInputs/1D_PML_centered/) with `final_time = 8.05` (~ParaView **cycle 80**, ~27 ns), while bulk fields and ψ are still O(1).

## Commands

```bash
cd testData/maxwellInputs/1D_PML_centered
./run_sign_matrix_t8.sh
python3 summarize_sign_matrix_t8.py
```

Artifacts: `sign_matrix_t8/S{0..6}.log`, `sign_matrix_t8/exports_S*`, `sign_matrix_t8/summary.txt`.

## Results at t = 8.05 (dt = 0.001, α = 0)

| Mode | `PML_SIGN_TEST` | Final `max\|psi\|` | `outer_Ey` | P0/P1 `\|Ey\|` @ end | `max\|Ey\|` probes | Gate (≤3) |
|------|-----------------|-------------------|------------|---------------------|-------------------|------------|
| default | 0 | **0.815** | 1.49 | 0.19 / 0.21 | 1.0 | **PASS** |
| S1 | 1 | 2.4×10¹⁴ | 8.9×10¹² | ~10⁹ | ~10⁹ | FAIL |
| S2 | 2 | 5.3×10⁶ | 2.7×10⁵ | 4066 / 5.1 | 7538 | FAIL |
| S3 | 3 | 4.1×10⁹ | 8.4×10⁷ | ~10⁶ | ~10⁶ | FAIL |
| S4 | 4 | 2.4×10¹⁴ | (same as S1) | (same as S1) | (same as S1) | FAIL |
| S5 | 5 | **0.815** | 1.49 | 0.19 / 0.21 | 1.0 | **PASS** |
| S6 | 6 | **0.815** | 1.49 | 0.19 / 0.21 | 1.0 | **PASS** |

**Criterion (interface window):** bounded `max|psi|` ~ O(1) and probe gate PASS.

## ψ timeline (default, selected Mult calls)

| Mult call | ~t_code | `max\|psi\|` | Notes |
|-----------|---------|--------------|-------|
| 27500 | 6.88 | 0.15 | Pre-peak |
| 30000 | 7.50 | **2.41** | Matches frame-75 audit |
| 31000 | 7.75 | 2.22 | |
| 31500 | 7.88 | 1.04 | |
| **32000** | **8.04** | **0.815** | End of short run |

Default run ends in a **post-peak ψ lull** (not monotonic growth through t = 8).

## Conclusions

1. **No sign mode improves** the interface-window state: S1/S4/S2/S3 are **strictly worse** even at t = 8.05; S5/S6 are **identical** to default (α = 0 → ψ mass inactive; global mult sign unchanged).
2. **Default at t = 8.05 is healthy** (`max|psi| < 1`, probes PASS) — confirms the trustworthy ParaView window (68–75 / ~23–25 ns) extends through ~cycle 80 in ψ, with a later lull before late blow-up (~t = 13+ in full runs).
3. **S2** (face driver +w) fails earlier than full t = 20 but with smaller ψ than S1/S3 — still unusable; face SBP opposition remains required.
4. Next work should stay on **formulation split / stiff integration** (sessions 15–18), not sign flips.

## Related

- [18-session-centered-vs-upwind-correlation.md](./18-session-centered-vs-upwind-correlation.md)
- [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md) (full t = 20 matrix)
