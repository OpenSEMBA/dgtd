# Cylindrical UPML `"R"` (Ji/Lu/Cai + face aux + Teixeira)

**Status: archive.** The `"R"` / CylUPML path is **removed** from the solver. JSON `active_axes: ["R"]` and `sigma_max` throw at parse. Cartesian `X`/`Y`/`Z` SC-PML is the live path — [`30-sc-pml-ade.md`](./30-sc-pml-ade.md). Do not treat this document as a live contract.

Preview-safe math: `$…$` / `$$…$$`.

---

## 1. Literature lock

| Source | Role |
|--------|------|
| Ji/Lu/Cai/Zhang JLT 2005 ([`refs/chen2005_dgtd_upml.pdf`](./refs/chen2005_dgtd_upml.pdf)) | UPML polarization ADE (TE+TM); unified phys/UPML state |
| Lu/Cai/Zhang JCP 2004 (cited [11]) | Parent Riemann / polarization discussion |
| Duru et al. ([`refs/duru2019_energy_stable_pml.pdf`](./refs/duru2019_energy_stable_pml.pdf)) | ADE auxiliaries **must** enter numerical fluxes |
| Teixeira/Chew | Cylindrical metric $\sigma_\theta=\Sigma/r$ |

**Rejected:** rotated Bagci SC-PML; Ji/Lu/Cai **volume-only** graft on vacuum Riemann; Kapidani maps; Bouquet centered UPML as primary.

---

## 2. Anti-reuse checklist

- No Bagci $a,b,c,d$; no `buildSCPMLOperators` wrapper
- No type/function name containing `SCPML` on the `"R"` path
- Own layout `CylUPMLLayout`, profiles `CylUPMLProfiles`, ops `cylUpmlVol_` / `cylUpmlFace_`
- Sharing allowed: FES, `MassIntegrator`, CSR merge helpers, Mult slot, MPI

---

## 3. Continuous equations (2D MVP, $\varepsilon=\mu=c=1$)

Principal frame $(\hat e_1,\hat e_2,\hat e_3)=(\hat r,\hat\theta,\hat z)$.

$$
\sigma_r=\sigma_r(\rho),\quad
\Sigma=\int_{r_{\mathrm{inner}}}^{r}\sigma_r\,dr',\quad
\sigma_1=\sigma_r,\quad
\sigma_2=\sigma_\theta=\Sigma/r,\quad
\sigma_3=0.
$$

Map Ji $(\sigma_x,\sigma_y)\mapsto(\sigma_1,\sigma_2)$. Combined TE+TM principal ADE (Ji eqs. (2)–(3)), then coefficient tensors to Cartesian via

$$
T_{xyz}=R\,\mathrm{diag}(T')\,R^{\mathsf{T}}
$$

**(geometry only — UPML $P,Q$ tensors, not Bagci).**

Aux packing in principal frame: $P=(Q^e_1,Q^e_2,P^m)$, $Q=(Q^m_1,Q^m_2,P^e)$.

Field damping principal diag: $(\sigma_2-\sigma_1,\sigma_1-\sigma_2,\sigma_1+\sigma_2)$.
$E\leftarrow P$ signs $(-1,-1,+1)$; $H\leftarrow Q$ signs $(+1,+1,-1)$.
Aux ODEs as in Ji (product drives + per-axis $\sigma$ damping on $Q^{e/m}$).

**TM-like (out-of-plane $E_3$ / in-plane $H$):** Ji (2) with $P^m,Q^m$.  
**TE-like (out-of-plane $H_3$ / in-plane $E$):** Ji (3) with $P^e,Q^e$.

State packing ($N=\mathrm{ndofs}$):

$$
[E_x E_y E_z\,|\,H_x H_y H_z\,|\,P_x P_y P_z\,|\,Q_x Q_y Q_z],
\quad n_{\mathrm{aux}}=6N.
$$

($P$/$Q$ store the packed TE+TM polarizations above after $R$-rotation of coefficient tensors.)

---

## 4. Discrete split (`GlobalEvolution`)

1. `globalOperator_->Mult` — curl + DG fluxes (**vacuum–PML faces ignored** when `"R"` active).
2. `cylUpmlFace_->AddMult` — Maxwell upwind **restored** on those faces **plus** polarization-dependent flux corrections (Duru/Lu–Cai).
3. `cylUpmlVol_->AddMult` — volume UPML ADE ($S^{(1)},S^{(2)}$ and TE analogues).

### Documented mods vs papers

1. Teixeira cylindrical $\sigma$ instead of slab $\sigma_x(x),\sigma_y(y)$.
2. MFEM assembled CSR instead of element-local DG.
3. Face aux: **additive** Duru-style $P/Q\to E/H$ injection on PML-side face DOFs of vacuum–PML faces (full Riemann replace deferred — interior faces lack attributes for `ignore_marker`). Scale $O(\alpha)$ upwind, not $O(\sigma)$. Vacuum Maxwell flux remains on those faces.

---

## 5. JSON

Required: `matches_vacuum`, `active_axes: ["R"]`, `sigma_max > 0`.  
Optional: `radial_center`, `grading_order` (default 0).  
Forbidden MVP: `kappa_max>1`, `alpha_max>0`.

---

## 6. Acceptance (2D Onion)

`2D_Onion_Dipole_R` / `_SMA`, `2D_Onion_TFSF_Scatter_R` / `_SMA`.

Gates: (0) no secular growth; (1) vacuum late $\|E\|$ ≪ SMA; (2) remanents may remain if decaying.

---

## 7. Session findings (2026-09-11) — wrap-up

### Shipped

Isolated `"R"` path was shipped then **removed** (not a Bagci wrapper). Historical map:

| Piece | Location |
|-------|----------|
| Props / JSON | `PMLProperties` — `cyl_upml`, `sigma_max`, `["R"]` alone |
| Profiles | `CylUPMLProfiles` — Teixeira $\sigma_r$, $\sigma_\theta=\Sigma/r$ |
| Layout | `CylUPMLLayout` — $n_{\mathrm{aux}}=6N$, packed $P,Q$ |
| Volume / face | `buildCylUPMLVolumeOperator`, `buildCylUPMLInterfaceFluxOperator` |
| Mult | `cylUpmlFace_` then `cylUpmlVol_` after `globalOperator_` |
| Cases | Onion R / R_SMA JSON already on schema; `docs/json-input-format.md` updated |

Cartesian `buildSCPMLOperators` / `scpmlOperator_` untouched.

### What failed earlier (do not resurrect)

1. **Rotated Bagci + Teixeira** — stable-ish but remanents / interface bounce.
2. **Ji volume-only on vacuum Riemann** — growth + sawteeth.
3. **First Cyl UPML draft** — wrong TE+TM blend (field damping $(\sigma_2+\sigma_3,\ldots)$ and wrong $P,Q$ signs/packing) → **exponential growth** on Dipole mid/shell. Face injection scaled by $O(\sigma)$ was also unstable at the interface with `grading_order=0`.

### Fix that made Dipole pass

Volume ADE locked to Ji JLT 2005 eqs. (2)–(3) with packing $P=(Q^e_1,Q^e_2,P^m)$, $Q=(Q^m_1,Q^m_2,P^e)$ and principal damping $(\sigma_2-\sigma_1,\sigma_1-\sigma_2,\sigma_1+\sigma_2)$. Face aux: additive $O(\alpha)$ on face DOFs (not $O(\sigma)$); vacuum Maxwell flux still on vac–PML faces (full ignore+Riemann deferred — interior faces lack `ignore_marker` attrs).

### Onion numbers (mpi-2, order 2)

| Case | Growth | Vacuum late vs SMA | Notes |
|------|--------|--------------------|-------|
| Dipole_R | None; mid/shell $\to\sim 10^{-6}$–$10^{-7}$ | $\sim 0.28\times$ SMA (**pass**) | Primary gate |
| TFSF_R | No exponential blow-up | $\sim 1.4$–$2\times$ SMA | Mild iface-out remanent / bounce; mid/shell stay small |

Exports: `Exports/mpi-2/2D_Onion_*_R*/PointProbes/`.

### Open for next session (not blockers for shutdown)

1. **Full vac–PML Riemann** — ignore those faces in `buildGlobalOperator`, restore Maxwell + aux flux in `cylUpmlFace_` (plan ideal; MVP is additive only).
2. **TFSF remanents** — improve late vacuum vs SMA; check iface-out half-moons in ParaView.
3. **Out of scope still:** 3D, CFS $\alpha$, $\kappa>1$, legacy `*PML_R*` RCS dirs, any Bagci+$R$ rotation.

### Success criteria check

- No `SCPML*` names on radial types — OK.
- No thin wrapper on `buildSCPMLOperators` — OK.
- Face path exists (not volume-only MVP) — OK (additive; full Riemann still open).
- Onion mid/shell do not grow — OK (Dipole); TFSF no secular blow-up.