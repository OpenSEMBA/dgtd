# AGENTS.md — OpenSEMBA/dgtd

Project guide for AI coding agents (Cursor, Claude Code, etc.). Behavioral rules live in [CLAUDE.md](./CLAUDE.md) and [.cursor/rules/](./.cursor/rules/); this file is **project-specific context**.

## What this repository is

**semba-dgtd** — time-domain Maxwell curl-equation solver using discontinuous Galerkin (DG) methods, built on a forked MFEM ([OpenSEMBA/mfem-geg](https://github.com/OpenSEMBA/mfem-geg)).

| Area | Location |
|------|----------|
| Main driver | `src/driver/` |
| Time evolution | `src/evolution/` — **`GlobalEvolution` is the active path** |
| DG operators | `src/components/DGOperatorFactory.h` |
| JSON → model | `src/driver/driver.cpp` |
| Cases | `testData/maxwellInputs/<case>/` |
| Tests | `test/` (gtest) |
| User docs | [docs/](./docs/) |

## Architecture (read before large changes)

### Evolution operators

| Operator | Status | Notes |
|----------|--------|-------|
| **`GlobalEvolution`** | **Active** | Assembled sparse DG operator + SGBC/TFSF; default in JSON |
| `MaxwellEvolution` | Deprecated | Do not add new features |
| `HesthavenEvolution` | Deprecated | Do not add new features |

Default JSON: `"evolution_operator": "global"` (or omit).

### GlobalEvolution operator split

1. **`globalOperator_`** — curl, DG fluxes, bulk conductivity (`DGOperatorFactory::buildGlobalOperator`).
2. **`Mult()` add-ons** — SGBC sub-solve + flux, TFSF source, (future) PML coupling.
3. **State vector** — `[Ex, Ey, Ez, Hx, Hy, Hz]` per DOF (`6 × ndofs`); PML will extend with auxiliary ψ (see PML docs).

### Units

Normalized internally: **c = ε = μ = 1**. SI only at I/O boundaries (`physicalConstants::*_SI`). See `src/components/Material.h`.

### Dimension-agnostic code

Always loop `X, Y, Z` and skip `d >= mesh.Dimension()`. Do not add 1D-only solver branches.

## Documentation map

| Topic | File |
|-------|------|
| **Build / run (human README)** | [README.md](./README.md) |
| **JSON input reference** | [docs/json-input-format.md](./docs/json-input-format.md) |
| **MOR → ParaView** | [docs/mor2paraview.md](./docs/mor2paraview.md) |
| **Volumetric CFS-CPML (in progress)** | [docs/pml/README.md](./docs/pml/README.md) |
| **MFEM/DGTD coding standards** | [.cursor/rules/02-mfem-dgtd-standards.mdc](./.cursor/rules/02-mfem-dgtd-standards.mdc) |
| **C++ guidance** | [.cursor/rules/03-cpp-expert-guidance.mdc](./.cursor/rules/03-cpp-expert-guidance.mdc) |

## Active work: volumetric CFS-CPML

Implementation plan and locked decisions: **`docs/pml/`** (start at [docs/pml/README.md](./docs/pml/README.md)).

Summary:

- **Formulation:** Gedney & Zhao ADE CFS-PML (DOI [10.1109/TAP.2009.2037765](https://doi.org/10.1109/TAP.2009.2037765)).
- **Integration:** `GlobalEvolution` only; attribute-tagged PML elements; extended ODE state; RK4 via `Mult()` first.
- **Remove:** surface `SBC_PML` shortcut (use volumetric PML only).
- **Branch:** `vol_pml` from `dev` (fresh implementation).

Before implementing PML, read `docs/pml/00-decisions-locked.md` and `docs/pml/04-implementation-plan.md`.

## Conventions for agents

1. **Minimal diffs** — match [CLAUDE.md](./CLAUDE.md); no drive-by refactors.
2. **JSON changes** — update [docs/json-input-format.md](./docs/json-input-format.md) and keep backward compatibility unless explicitly removing a feature.
3. **Commits** — only when the user asks.
4. **Tests** — run relevant gtests after substantive changes; PML acceptance is manual via `testData/maxwellInputs/` (−40 dB reflection target).
5. **MFEM** — must use `external/mfem-geg` fork, not upstream MFEM.

## Quick build (MPI release)

```sh
export METIS_DIR=...   # see README.md
export HYPRE_DIR=...
cmake --preset gnu-release-mpi
cmake --build --preset build-gnu-release-mpi --parallel
```

Binary typically: `build/gnu-release-mpi/bin/opensemba_dgtd` (confirm in your preset output directory).

## Cursor-specific notes

Cursor loads:

- **Always-applied:** `.cursor/rules/01-core-behavior.mdc`
- **Glob-triggered:** `02-mfem-dgtd-standards.mdc`, `03-cpp-expert-guidance.mdc` on `src/**`, `test/**`
- **This file:** `AGENTS.md` — read when onboarding to the repo or starting multi-file features

There is no separate Cursor equivalent to `CLAUDE.md` beyond rules + `AGENTS.md`; keep `AGENTS.md` updated when architecture or doc locations change.
