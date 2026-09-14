# Radial DGTD–PML — literature vs live path

**Status: archive.** The `"R"` / CylUPML path is **removed** from the solver. JSON `active_axes: ["R"]` throws at parse. Notes below describe the former implementation — [`33-radial-sc-pml-ade.md`](./33-radial-sc-pml-ade.md).

**Forbidden to resurrect:** rotated Bagci + vacuum faces; volume-only Lu/Cai graft.

| File | Paper | Relevance |
|------|-------|-----------|
| `chen2005_dgtd_upml.pdf` | Ji/Lu/Cai/Zhang JLT 2005 | UPML ADE TE+TM; DGTD |
| `duru2019_energy_stable_pml.pdf` | Duru/Gabriel/Kreiss | Aux must enter face fluxes |
| `feizhou2008_dg_pml.pdf` | Bouquet et al. | DG+UPML (centered; reference only) |
| `brivio2020_matrixfree_pml.pdf` | Kapidani/Schöberl | Stretch maps — parked |

## Cylindrical metric (Teixeira/Chew)

$$
s_r = 1 + \frac{\sigma_r}{j\omega},\quad
\tilde r = r_{\mathrm{inner}} + \int_{r_{\mathrm{inner}}}^{r} s_r\,dr',\quad
s_\theta = \frac{\tilde r}{r},\quad
\sigma_\theta = \frac{\Sigma}{r}.
$$
