# Radial DGTD–PML — literature vs live path

**Live `"R"` path (2026-09):** **cylindrical SC-PML** in 2D — $\sigma_r(\rho)$ graded, $\sigma_\theta=\Sigma/r$ with $\Sigma=\int_{r_{\mathrm{inner}}}^{r}\sigma_r\,dr'$, $\sigma_z=0$, $\kappa\equiv 1$; Bagci $a,b,c,d$ in $(\hat r,\hat\theta,\hat z)$; $T=R\,\mathrm{diag}\,R^T$; volume `Mass` ADE only; vacuum upwind flux unchanged (`matches_vacuum`).

**Superseded:** locally uniaxial quasi-PML with $\sigma=(\sigma_r,0,0)$ (worse oblique / stronger $E_r$ issues).

**Deferred:** 3D spherical $(s_\theta,s_\phi)$; CFS $\alpha$; Duru face penalties on $P$.

Note: with $\kappa\equiv 1$ and $\sigma_z=0$, Bagci still has a DC null space on $E_r$ and $E_\theta$ ($\sigma_v\sigma_w=0$); cylindrical $\sigma_\theta$ helps propagating waves but does not remove that algebraic null space. See [`30-sc-pml-ade.md`](./30-sc-pml-ade.md) §4.4.

| File | Paper | Relevance |
|------|-------|-----------|
| `chen2005_dgtd_upml.pdf` | Ji/Lu/Cai/Zhang JLT 2005 | DGTD + **Cartesian UPML** unified state; flux at phys/PML |
| `feizhou2008_dg_pml.pdf` | Bouquet/Dedeban/Piperno HAL 2008 | DG + UPML; centered flux; doubled unknowns in PML |
| `brivio2020_matrixfree_pml.pdf` | Kapidani/Schöberl arXiv:2002.08733 | Matrix-free DG + **complex stretch** via covariant maps |
| `duru2019_energy_stable_pml.pdf` | Duru/Gabriel/Kreiss arXiv:1802.06388 | Energy-stable DGSEM PML; **fluxes must extend to ADE auxiliaries** |

Wu et al. BoR-DGTD (2023) is paywalled — not downloaded.

---

## Architecture (unchanged skeleton)

Chen-style split: hyperbolic curl/flux stays vacuum (`globalOperator_`); PML damping/ADE is volume-only (`scpmlOperator_`). Core–PML faces keep ordinary upwind with `matches_vacuum`.

## Cylindrical metric (Teixeira/Chew)


$$
s_r = 1 + \frac{\sigma_r}{j\omega},\quad
\tilde r = r_{\mathrm{inner}} + \int_{r_{\mathrm{inner}}}^{r} s_r\,dr',\quad
s_\theta = \frac{\tilde r}{r} = 1 + \frac{\Sigma}{j\omega\, r},\quad
\sigma_\theta = \frac{\Sigma}{r}.
$$


Power-law grading $\sigma_r=\sigma_{\max}(\rho/L)^m$ ⇒ $\Sigma=\sigma_{\max} L/(m+1)\,(\rho/L)^{m+1}$.

## Historical notes

- CFS-lite ($d=\sigma+\alpha$ only) grew fields / hung — abandoned.
- Probe remanents under quasi-UPML tracked $E_r$ (rotation OK; formulation incomplete).
- Duru auxiliary-face fluxes remain a possible discrete-stability follow-up if late junk survives cylindrical $\sigma_\theta$.
