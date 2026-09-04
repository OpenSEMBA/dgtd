# 1D_PML_L3

Same as [`1D_PML_L2`](../1D_PML_L2/) with PML thickness **L = 3** on \(x\in[2,5]\)
(30 elements), buffer to \(x=5.2\).

### 2026-09-04 parameter sweep (probe 0, inc [3,7.2], ref [8.5,11.5])

| L | m | R_target | R_dB(f_peak) | gate |
|---|---|----------|--------------|------|
| 1 (`1D_PML_DFT`) | 3 | 1e-4 | −38.5 | FAIL |
| 1 | 3 | 1e-6 | −31.9 | FAIL |
| 1 | 3 | 1e-8 | −28.0 | FAIL |
| 1 | 4 | 1e-4 | −34.7 | FAIL |
| 1 | 4 | 1e-6 | −28.8 | FAIL |
| 1 | 4 | 1e-8 | −25.4 | FAIL |
| **2** | 3–4 | 1e-4…1e-8 | **−148 … −161** | **PASS** |
| **3** | 3–4 | 1e-4…1e-8 | **−178 … −191** | **PASS** |

Takeaways:

1. **Thickness dominates.** Doubling L (1→2) crosses −40 dB by a large margin.
2. On the thin (L=1) mesh, **more conductivity (smaller R_target) made reflection worse** —
   discrete interface error, not σ budget.
3. Milder grading (`m=3`) was slightly better than `m=4` at L=1, irrelevant once L≥2.
4. All L≥2 runs stayed late-time stable (`max|ψ|` at t=20 ~ 1e−14).
