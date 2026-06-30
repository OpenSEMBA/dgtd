# Time integration: RK4, Mult, ImplicitSolve

## MFEM contract

`mfem::TimeDependentOperator` provides:

| Method | Used by | Meaning |
|--------|---------|---------|
| `Mult(x, y)` | Explicit solvers (RK4) | Compute **y = F(x)** = time derivative **dx/dt** |
| `ImplicitSolve(dt, x, k)` | Implicit / SDIRK solvers | Solve implicit stage equations involving **I - dt*A** |

dgtd `GlobalEvolution` implements both; default user setting is **RK4**.

---

## Milestone A: explicit RK4 only

### Locked decision

- All CFS-CPML dynamics (E, H, **ψ**) integrated **explicitly** via **`Mult()`**.
- Same RK4 stages for vacuum and PML regions.
- **No** regional IMEX (explicit vacuum / implicit PML).
- **No** `ImplicitSolve` changes required for Milestone A pass.

### Extended state in RK4

ODE vector:

```text
x = [ E_H (6N) ; ψ (n_aux) ]
```

Each RK4 stage calls `Mult(x_stage, k)` where **k** has same size as **x**:

```text
k = [ dE/dt, dH/dt ; dψ/dt ]
```

ψ **must** participate in every stage — **not** updated once per step like SGBC.

### SGBC interaction (unchanged)

`Solver::step()`:

```cpp
commitSGBCCheckpoint(...);
odeSolver_->Step(fields_.allDOFs(), ...);  // RK4 calls Mult multiple times
finalizeSGBCStep(...);
```

SGBC checkpoint wraps the **whole** RK4 step. PML ψ is **inside** `allDOFs()` and advances with each internal `Mult()` — **no** SGBC-style checkpoint for ψ.

---

## Why not regional IMEX for v1

User idea: explicit vacuum + implicit PML.

| Aspect | Assessment |
|--------|------------|
| Physics | Valid split if curl explicit, stiff α ψ implicit |
| MFEM | No built-in regional IMEX; custom solver or manual splitting |
| dgtd today | `ImplicitSolve` implicitizes **entire** `globalOperator_`, not PML block only |
| Effort | High — new time integration layer |
| Need | Try explicit first; user accepted RK4 for first test |

**Defer** regional IMEX to possible future Milestone C.

---

## Current `ImplicitSolve` behavior (reference)

`GlobalEvolution::ImplicitSolve(dt, x, k)`:

1. Asserts `x.Size() == 6 * ndofs` (**will break** when extended — not used in Milestone A).
2. Computes `rhs = globalOperator_ * x` (+ TFSF).
3. Solves `(I - dt*A) k = rhs` via dense LU (small) or GMRES.

**Entire** semidiscrete operator treated implicitly — includes curl and bulk σ, **not** PML aux today.

---

## Milestone B: SDIRK option (when explicit insufficient)

### Trigger

- Late-time instability with `alpha_max > 0` still insufficient on RK4
- Need larger dt than explicit stability allows for ψ equations

### Approach (simpler than regional IMEX)

When model `hasPML()`:

- Auto-select **`SDIRK23`** in `Solver::assignODESolver()` (comment already mentions PML in `Solver.cpp`).
- Extend **`ImplicitSolve`** to size `6N + n_aux`.
- Split operator:

```text
Explicit part (Mult): curl flux + ADE ψ coupling as currently written
Implicit part (ImplicitSolve): treat stiff terms — at minimum α ψ decay block
```

**Full implicit monolithic** `(I - dt*A_full)` including curl is also possible but expensive; **stiffly stable split** targeting ψ and σ coupling is preferred.

### MFEM solvers available (`SolverOptions.h`)

| Enum | Type |
|------|------|
| RK4 | Explicit |
| BackwardEuler, Trapezoidal, ImplicitMidpoint | Implicit |
| SDIRK33, SDIRK23, SDIRK34 | Implicit RK |

User locked **RK4 first**; SDIRK is Milestone B.

---

## Stability heuristics (normalized units)

### CFL (vacuum / curl)

Existing `estimateTimeStep()` in `Solver.cpp` — unchanged for Milestone A.

### PML auxiliary stiffness

ADE terms **∂ψ/∂t = −α ψ + …** suggest explicit stability roughly **α · Δt ≲ const** (order 1–2 depending on RK order).

Monitor:

- If `alpha_max` large relative to `1/dt`, expect explicit issues
- Mitigation: increase grading depth, reduce `alpha_max`, reduce `dt`, or enable SDIRK (B)

### CFS purpose

**α > 0** improves late-time behavior of PML **medium**; distinct from time-stepper choice.

---

## Explicit vs bulk conductivity

| Mechanism | Integration today | PML? |
|-----------|-------------------|------|
| Bulk σ (Material) | In `globalOperator_`, explicit under RK4 | **No** — forbidden on PML tags |
| PML stretch σ | ADE via `applyPMLCoupling` | **Yes** |

Do not conflate the two when extending `ImplicitSolve`.

---

## Implementation notes for `Mult()` ordering

Recommended order (consistent with today):

```text
1. SGBC sub-solve (if any)
2. globalOperator_->Mult  → curl + κ-mass + bulk σ on E
3. applyPMLCoupling       → ADE ψ and corrections
4. TFSF source
5. SGBC flux injection
```

Document final order in code comment when implemented. TFSF before/after PML may matter for TF regions adjacent to PML — user stated no overlap.

---

## Summary table

| Phase | Integrator | PML ψ | ImplicitSolve |
|-------|------------|-------|---------------|
| Milestone A | RK4 | Explicit in Mult | Unused / broken size OK if RK4 only |
| Milestone B | SDIRK23 optional | Explicit or split implicit | Extended + PML blocks |
| Future | Regional IMEX | TBD | Custom |

---

## Code touch list for Milestone B (preview only)

- `GlobalEvolution::ImplicitSolve` — extended sizes, include `PMLOperator_` in A matrix
- `Solver::assignODESolver` — auto SDIRK when `model.hasPML()`
- `SolverOptions` / JSON — optional `ode_type` override documented

Not in scope until Milestone A acceptance.
