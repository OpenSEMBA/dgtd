# Conversation context snapshot

This file captures planning decisions from the Cursor design conversation so future sessions need not reload chat history.

**Date context:** Planning completed before first implementation Agent pass.  
**User:** Alejandro (OpenSEMBA/dgtd).  
**Advisor:** Salvador González García (UGR / GEG) — ADE CFS-PML via Gedney & Zhao paper + Taflove CPML intuition.

---

## User constraints stated explicitly

1. **GlobalEvolution only** — Maxwell/Hesthaven deprecated for new work.
2. **Dimension-agnostic** — cannot add 1D-only code path; 1D mesh is test config only.
3. **Gmsh attributes** define PML regions; thickness from mesh, not JSON.
4. **No VolumetricRegionSubMesher** required for v1.
5. **Remove SBC_PML** and old vol PML scraps; recreate **`vol_pml`** branch from **`dev`**.
6. **CPU, serial, single mesh** for Milestone A.
7. **SGBC/TFSF** coexist but never on PML interfaces in practice.
8. **No automated tests** — manual `testData/maxwellInputs` cases.
9. **Probes in vacuum**; user DFT offline; **−40 dB** acceptance.
10. **Normalized units** throughout.
11. **RK4 explicit first**; ImplicitSolve later is fine.
12. **ψ not exported** to Paraview.
13. **matches_vacuum** — no bulk conductivity near PML.
14. **Outer boundary** — implementer chooses SMA over exploratory PEC outers.

---

## Key clarifications resolved

### CPML vs CFS-PML vs UPML

- **CFS-PML** = stretch with α in denominator.
- **CPML** = time-domain via convolution or equivalent **ADE**.
- **Implementation name:** **ADE CFS-CPML** (Gedney & Zhao 2010).
- Not contradictory to recommend auxiliary ODEs — that **is** CPML in DGTD form.

### IMEX vacuum explicit / PML implicit

- Valid physics, **deferred** — overkill for v1.
- Full **SDIRK on extended system** is simpler fallback before regional IMEX.

### `pml_thickness` JSON confusion

- User always meshed PML volume in Gmsh with attribute tags.
- JSON controls grading law and max parameters, **not** geometric thickness.

---

## Agent instructions for next session

1. Read `docs/pml/README.md` first.
2. Follow `04-implementation-plan.md` step order (A.1 → A.6).
3. Respect `00-decisions-locked.md` — do not re-debate.
4. Transcribe Maxwell ADE equations from DOI **10.1109/TAP.2009.2037765** when coding A.5.
5. Stop after each step's verification criteria; ask user before long physics debug loops.
6. Do not commit unless user asks.

---

## Files created in docs/pml/

| File |
|------|
| README.md |
| 00-decisions-locked.md |
| 01-physics-and-formulation.md |
| 02-codebase-architecture.md |
| 03-json-schema-and-mesh.md |
| 04-implementation-plan.md |
| 05-verification.md |
| 06-cleanup-inventory.md |
| 07-time-integration.md |
| 08-conversation-snapshot.md (this file) |

---

## Latest session (2026-05-27)

Implementation log: **[09-session-record-2026-05-27.md](./09-session-record-2026-05-27.md)**. Completed Milestone A.5 (factory `PMLOperator_`, volume+face ψ driver, component-indexed vector ψ). A.6 partially verified on `1D_PML/` — resume with user DFT and late-time SMA-edge tuning (`alpha_max`, SDIRK).

---

## Documentation layout (updated)

- **README.md** — repository overview and build/run only.
- **AGENTS.md** — AI agent project context (architecture, doc map, PML pointer).
- **docs/json-input-format.md** — full JSON reference (moved from README).
- **docs/mor2paraview.md** — MOR post-processing (moved from README).
- **docs/pml/** — CFS-CPML design memory.
