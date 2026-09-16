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
2. **`Mult()` add-ons** — SGBC sub-solve + flux; Cartesian SC-PML ADE (`scpmlCurlDelta_` then `scpmlOperator_`); TFSF source.
3. **State vector** — `[Ex, Ey, Ez, Hx, Hy, Hz]` per DOF (`6 × ndofs`); with Cartesian PML, plus Bagci $P$ auxiliaries (`SCPMLLayout`, $n_{\mathrm{aux}}=6N$).

### Units

Normalized internally: **c = ε = μ = 1**. SI only at I/O boundaries (`physicalConstants::*_SI`). See `src/components/Material.h`.

### Dimension-agnostic code

Always loop `X, Y, Z` and skip `d >= mesh.Dimension()`. Do not add 1D-only solver branches.

## Documentation map

| Topic | File |
|-------|------|
| **Build / run (human README)** | [README.md](./README.md) |
| **JSON input reference** | [docs/json-input-format.md](./docs/json-input-format.md) (Maxwell cases + offline `opensemba_rcs`) |
| **MOR → ParaView** | [docs/mor2paraview.md](./docs/mor2paraview.md) |
| **Volumetric PML (SC-PML ADE)** | [docs/pml/README.md](./docs/pml/README.md) |
| **MFEM/DGTD coding standards** | [.cursor/rules/02-mfem-dgtd-standards.mdc](./.cursor/rules/02-mfem-dgtd-standards.mdc) |
| **C++ guidance** | [.cursor/rules/03-cpp-expert-guidance.mdc](./.cursor/rules/03-cpp-expert-guidance.mdc) |

## Active work: volumetric PML (Cartesian SC-PML)

- **Cartesian** `X`/`Y`/`Z`: Bagci/Chen SC-PML — [`docs/pml/30-sc-pml-ade.md`](./docs/pml/30-sc-pml-ade.md). JSON `"R"` is rejected at parse.
- Archived radial notes: [`33`](./docs/pml/33-radial-sc-pml-ade.md), [`32`](./docs/pml/32-radial-dgtd-literature.md). Survey / status: [`29`](./docs/pml/29-dgtd-pml-approaches.md), [`31`](./docs/pml/31-sc-pml-status.md).
- **Integration:** `GlobalEvolution` only. Cartesian volume ADE after `globalOperator_`. Do **not** reintroduce Gedney $\psi\sim D(F)$, Bagci+$R$ rotation, or the deleted CylUPML path.

## Conventions for agents

1. **Minimal diffs** — match [CLAUDE.md](./CLAUDE.md); no drive-by refactors.
2. **JSON changes** — update [docs/json-input-format.md](./docs/json-input-format.md) and keep backward compatibility unless explicitly removing a feature.
3. **Commits** — only when the user asks.
4. **Tests** — run relevant gtests after substantive changes; PML acceptance is manual via `testData/maxwellInputs/` (−40 dB reflection target).
5. **MFEM** — must use `external/mfem-geg` fork, not upstream MFEM.
6. **NEVER clean builds/MFEM** — do not `ninja clean`, wipe `build/`, or delete MFEM/CUDA `.o` / `libmfem.a`. Incremental rebuild only. See [.cursor/rules/04-never-clean-build.mdc](./.cursor/rules/04-never-clean-build.mdc).

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

- **Always-applied:** `.cursor/rules/01-core-behavior.mdc`, `04-never-clean-build.mdc`, `05-markdown-math-preview.mdc`, `06-frozen-contracts.mdc`
- **Glob-triggered:** `02-mfem-dgtd-standards.mdc` and `03-cpp-expert-guidance.mdc` on `src/` and `test/` C++/CUDA; `07-deprecated-evolution.mdc` on deprecated evolution + Hesthaven tests; `08-json-tag-contract.mdc` on driver/JSON; `09-verification-policy.mdc` on `test/`
- **Frozen contracts:** [`.cursor/contracts/frozen-contracts.md`](.cursor/contracts/frozen-contracts.md)
- **Specialist agents:** [`.cursor/agents/`](.cursor/agents/) — Team Lead plus Formulation, FaceFlux, OperatorAssembly, Evolution, Execution, MFEM, Verification; separate PML and SGBC; scoped Configuration, TFSF/Sources, Probes/RCS/Export; literature Bookkeeper and Lore Interpreter (advisory)
- **Commands:** `/dgtd-intake`, `/dgtd-consensus`, `/dgtd-literature`, `/dgtd-review`
- **This file:** `AGENTS.md` — read when onboarding or starting multi-file features

Do **not** migrate `.opencode` into `.cursor`. `.opencode/instructions/mfem-dgtd.md` is incomplete evidence (it mentions a nonexistent `core/schema` / `.dgtd.json`). Live JSON schema is [docs/json-input-format.md](docs/json-input-format.md).

There is no separate Cursor equivalent to `CLAUDE.md` beyond rules + `AGENTS.md`; keep `AGENTS.md` updated when architecture or doc locations change.

### Multi-agent protocol (short)

Classify blast radius (LOCAL / PAIR / SPATIAL-TRIAD / FULL). Assign the minimum specialists. No majority vote on signs, flux, layout, or device residency. No coupled numerical implementation without frozen-contract confirmation **and explicit user approval**. After implementation, Verification plus owning specialists review. A paper does not override repository evidence, tests, MFEM, or a user decision. Do not extend `MaxwellEvolution` or `HesthavenEvolution` (`HesthavenEvolution` is a read-only operator oracle).
