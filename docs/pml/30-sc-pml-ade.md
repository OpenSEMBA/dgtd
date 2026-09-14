# SC-PML ADE (Bagci/Chen) — formulation ↔ code

**Status: active** (replaces CuDG3D classical $J$/$M$ and paused Gedney $\psi\sim D(F)$).

This note is the **paper-style** description of what dgtd actually does: continuous equations, how they become a sparse RHS operator, where rotation enters, and which files/functions implement each step. Read it before changing PML assembly or `Mult`.

**CFS $\alpha$:** not wired. JSON must keep `alpha_max: 0` (parser rejects otherwise).

> **Preview note:** display math uses `$$ ... $$` with no line-leading `-` inside blocks, so Cursor/VS Code KaTeX preview is not broken by Markdown list parsing.

---

## 1. Continuous equations (principal / Cartesian axes)

Normalized units $\varepsilon=\mu=c=1$. For each field component $u\in\{x,y,z\}$ with cyclic partners $(u,v,w)$:

**(1)**

$$
\begin{aligned}
\partial_t (a_{uu} E_u)
  &= (\nabla\times H)_u - b_{uu}\, E_u - c_{uu}\,(P_E)_u, \\
\partial_t (a_{uu} H_u)
  &= -(\nabla\times E)_u - b_{uu}\, H_u - c_{uu}\,(P_H)_u, \\
\partial_t (P_E)_u
  &= \kappa_u^{-1}\, E_u - d_{uu}\,(P_E)_u, \\
\partial_t (P_H)_u
  &= \kappa_u^{-1}\, H_u - d_{uu}\,(P_H)_u.
\end{aligned}
$$

Diagonal stretch tensors (inactive axis: $\sigma=0$, $\kappa=1$):

**(2)**

$$
\begin{aligned}
a_{uu} &= \frac{\kappa_v\kappa_w}{\kappa_u}, \\
b_{uu} &= \frac{\sigma_v\kappa_w + \sigma_w\kappa_v - a_{uu}\sigma_u}{\kappa_u}, \\
c_{uu} &= \sigma_v\sigma_w - b_{uu}\sigma_u, \\
d_{uu} &= \frac{\sigma_u}{\kappa_u}.
\end{aligned}
$$

At $\kappa\equiv 1$, $a_{uu}=1$ and (1)–(2) reduce to the CuDG3D $J=\sigma^2 P$ classical ADE.

### 1.1 What “stretch in $X$” actually damps

Stretch profiles $(\sigma_u,\kappa_u)$ live on **coordinate axes**, not on “the field named $E_u$” alone. For a **uniaxial** slab with only $\sigma_x>0$ and $\sigma_y=\sigma_z=0$, $\kappa\equiv 1$:

| Component $u$ | $b_{uu}$ | $c_{uu}$ | $d_{uu}$ | Field ODE (aside from curl) |
|---------------|----------|----------|----------|-----------------------------|
| $x$ | $-\sigma_x$ | $\sigma_x^2$ | $\sigma_x$ | $\partial_t E_x \mathrel{+}= +\sigma_x E_x - \sigma_x^2 (P_E)_x$ |
| $y$ | $+\sigma_x$ | $0$ | $0$ | $\partial_t E_y \mathrel{+}= -\sigma_x E_y$ |
| $z$ | $+\sigma_x$ | $0$ | $0$ | $\partial_t E_z \mathrel{+}= -\sigma_x E_z$ |

So an $X$-stretch:

- puts the **ADE memory** $(P_E)_x,(P_H)_x$ on the **normal** components $E_x,H_x$;
- applies **instantaneous damping** $-\sigma_x$ to the **tangential** components $E_y,E_z$ (and $H_y,H_z$).

It does **not** mean “kill $E_x$ only.” Tangential fields are damped without auxiliaries; the normal field is handled by the $(E,P)$ pair. Same pattern for a $Y$- or $Z$-only slab by cycling indices.

**Multi-axis corners** (`active_axes: ["X","Y"]`, …): several $\sigma_u$ nonzero → $c_{uu}$ couples auxiliaries biaxially (true corner stretch, not a hacked sum of 1D slabs).

---

## 2. Discrete state and operator split

### 2.1 State vector

When any PML region exists, the RK4 unknown is the **extended** vector
(`SCPMLLayout`, [`src/components/SCPMLLayout.h`](../../src/components/SCPMLLayout.h)):

$$
\mathbf{v}
=
\bigl[
  E_x,\,E_y,\,E_z,\,
  H_x,\,H_y,\,H_z
  \;\big|\;
  (P_E)_x,\,(P_E)_y,\,(P_E)_z,\,
  (P_H)_x,\,(P_H)_y,\,(P_H)_z
\bigr]^{\mathsf{T}}
\in\mathbb{R}^{12N}.
$$

with $N=$ local DOFs. Offsets: fields $0\ldots 6N-1$; $P_E$ at $6N+cN$; $P_H$ at $9N+cN$ (`pEOffset` / `pHOffset`).

If there is no PML, $\mathbf{v}$ is only the $6N$ fields and `scpmlOperator_` is null.

### 2.2 RHS as two operators (not “mystery adds”)

Time marching integrates $\dot{\mathbf{v}} = \mathbf{f}(\mathbf{v},t)$. Inside `GlobalEvolution::Mult` ([`src/evolution/GlobalEvolution.cpp`](../../src/evolution/GlobalEvolution.cpp)):

**(3)**

$$
\mathbf{f}
=
A_{\mathrm{EH}}\,
\mathbf{v}_{\mathrm{EH}}^{\mathrm{(nbr)}}
{+}
A_{\mathrm{PML}}\,
\mathbf{v}
{+}
\text{(TFSF / SGBC sources)}.
$$

| Symbol | Role | Code |
|--------|------|------|
| $A_{\mathrm{EH}}$ | DG curl + interior/boundary fluxes, **vacuum** masses; width includes face-neighbour ghosts | `globalOperator_->Mult(multWorkVec_, out_fields)` into the first $6N$ slots of `out` |
| $\mathbf{v}_{\mathrm{EH}}^{\mathrm{(nbr)}}$ | Local EH DOFs + exchanged face-neighbour EH | built into `multWorkVec_` before the Mult |
| $A_{\mathrm{PML}}$ | Volume ADE only (damping $b,c$ + $P$ ODEs); **no** face fluxes | `scpmlOperator_->AddMult(in, out)` |
| optional $\Delta_u$ | If $\kappa_{\max}>1$: replace unit-mass curl by $M_a^{-1}$ weak form on PML DOFs | `scpmlCurlDelta_[u]` after the EH Mult |

So PML is **not** folded into `globalOperator_`. It is a second sparse matrix acting on the **full** extended state and **added** into the same RHS `out`:

```text
out ← 0
out[0:6N)  ← A_EH * v_EH_nbr          // curl/flux
out[0:6N)  += Δ * out[0:6N)           // optional κ>1 only
out        += A_PML * v               // ADE (fields + P)
```

Equivalent statement: $\dot{\mathbf{v}} = A_{\mathrm{total}}\mathbf{v}+\mathbf{s}$ with
$A_{\mathrm{total}} = \mathrm{blkdiag}(A_{\mathrm{EH}},0_{P}) + A_{\mathrm{PML}}$
(plus the $\kappa>1$ curl post-process).

---

## 3. How $A_{\mathrm{PML}}$ is built

Factory: `DGOperatorFactory::buildSCPMLOperators` in
[`src/components/DGOperatorFactory.h`](../../src/components/DGOperatorFactory.h).

### 3.1 Weak form ($\kappa\equiv 1$)

On PML-marked elements, unit inverse mass $M^{-1}$ times graded masses:

**(4)**

$$
\begin{aligned}
\dot{\mathbf{E}}
  &\mathrel{+}=
  {-}\,M^{-1} M(b)\,\mathbf{E}
  {-}\,M^{-1} M(c)\,P_E, \\
\dot{P}_E
  &\mathrel{+}=
  {+}\,M^{-1} M(\kappa^{-1})\,\mathbf{E}
  {-}\,M^{-1} M(d)\,P_E.
\end{aligned}
$$

and the same pattern for $H,P_H$. Here $M(\phi)$ is the DG mass with scalar (or, for radial, tensor-entry) coefficient $\phi$ assembled **only on PML volume attributes**.

### 3.2 Cartesian path (diagonal $b_{uu}$, …)

For each direction $u\in\{X,Y,Z\}$:

1. `SCPMLTensorCoefficient` evaluates (2) at each QP from
   `PMLProfileData::evaluateAtTransform` ($\sigma,\kappa$ on that axis).
2. Marked `MassIntegrator` → $M(b_u)$, $M(c_u)$, $M(d_u)$, $M(1/\kappa_u)$.
3. Left-multiply by `MInv[E]` / `MInv[H]` (`buildByMult`).
4. Place blocks into the big CSR with **signs matching (1)/(4)**:

| Continuous term | Block $A_{\mathrm{PML}}$ | Scale in `collectBlockPlacement` |
|-----------------|--------------------------|----------------------------------|
| $-b_{uu} E_u$ | row $E_u$, col $E_u$ | $-1$ |
| $-c_{uu}(P_E)_u$ | row $E_u$, col $(P_E)_u$ | $-1$ |
| $\kappa_u^{-1} E_u$ | row $(P_E)_u$, col $E_u$ | $+1$ |
| $-d_{uu}(P_E)_u$ | row $(P_E)_u$, col $(P_E)_u$ | $-1$ |

(and identical for $H$). **No off-diagonal Cartesian field coupling** when axes are mesh-aligned: $A_{\mathrm{PML}}$ is $3$ independent $2\times 2$ (field,aux) pairs per polarization, times $E$ and $H$.

### 3.3 $\kappa>1$ (Cartesian only today)

- Damping uses $M_a^{-1}M(b)$, $M_a^{-1}M(c)$ instead of $M^{-1}M(\cdot)$.
- After $A_{\mathrm{EH}}$ Mult, each field component slice is corrected by
  $\Delta_u = M_a^{-1}M - M^{-1}M$ on PML DOFs so the curl piece matches LHS $a_{uu}$.

---

## 4. Radial / cylindrical path (`active_axes: ["R"]`) — removed

Cylindrical `"R"` (Teixeira metrics + rotated Bagci, and later Lu/Cai grafts) is **removed**. JSON `["R"]`/`["r"]` throws at parse. Archive: [`33-radial-sc-pml-ade.md`](./33-radial-sc-pml-ade.md).

Optional `stretch_mode: "radial"` still means **radial depth grading** for Cartesian `X`/`Y`/`Z` ADE stacks only (isotropic $\rho$; not $\sigma_\theta=\Sigma/r$).

---

## 5. End-to-end map (checklist)
| Step | Math | Code |
|------|------|------|
| Parse PML JSON | regions, axes, grading | `parsePMLMaterialBlock` — [`PMLProperties.cpp`](../../src/components/PMLProperties.cpp) |
| $\sigma,\kappa$ at QPs | depth $\rho$, power law | `PMLProfileData` — [`PMLProfiles.cpp`](../../src/components/PMLProfiles.cpp) |
| Layout $P_E,P_H$ | $n_{\mathrm{aux}}=6N$ | `SCPMLLayout` |
| Build $A_{\mathrm{PML}}$ | (4), CSR merge | `buildSCPMLOperators` |
| Build $A_{\mathrm{EH}}$ | DG Maxwell | `buildGlobalOperator` (unchanged by PML) |
| RHS | (3) | `GlobalEvolution::Mult` |
| Time integrate | RK4 on $\mathbf{v}$ | `odeSolver_->Step` |

---

## 6. JSON (live)

No `pml_formulation` switch. On `"type": "PML"`:

| Field | Default | Notes |
|-------|---------|-------|
| `kappa_max` | `1` | ≥1; curl $a$-rescale for Cartesian |
| `alpha_max` | `0` | Must stay 0 (CFS deferred; parser rejects `>0`) |
| `active_axes` | required | `"X"`/`"Y"`/`"Z"` only (`["R"]` throws at parse — archive [`33`](./33-radial-sc-pml-ade.md)) |
| `grading_order`, `target_reflection`, `stretch_mode`, `radial_center` | as before | `stretch_mode: "radial"` = radial depth for Cartesian axes |

Cartesian only. Onion `"R"` cases on disk do not run.

---

## 7. Acceptance (practical)

| Case | Gate |
|------|------|
| `1D_PML` | DFT ≤ −40 dB ($\kappa=1$) |
| `2D_PML_X_slab` | Late-stable |
| Dipole Cartesian PML vs SMA | Vacuum bounce ≪ SMA |

---

## 8. References

- Chen et al., arXiv:2006.02551 — SC-PML ADE  
- Teixeira & Chew — cylindrical/spherical metric $s_\theta=\tilde r/r$  
- Survey: [`29-dgtd-pml-approaches.md`](./29-dgtd-pml-approaches.md)  
- Radial notes: [`32-radial-dgtd-literature.md`](./32-radial-dgtd-literature.md)  
- Status / history: [`31-sc-pml-status.md`](./31-sc-pml-status.md)  
- Gedney archive: [`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md)
