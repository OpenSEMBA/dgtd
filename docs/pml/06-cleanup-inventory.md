# Cleanup inventory: remove before / during Milestone A

## Policy

**Start from zero** for volumetric PML. Remove surface PML shortcut and dead experimental hooks. Do **not** merge old `vol_pml` branch code.

---

## Remove: `SBC_PML` (surface absorber shortcut)

### `src/driver/driver.cpp`

| Symbol / region | Action |
|-----------------|--------|
| `isSGBCBoundaryType()` | Remove `\|\| boundary_type == "SBC_PML"` |
| `buildAutoPMLSGBCLayers` lambda (~1431–1478) | **Delete entire function** |
| Boundary loop `if (boundary_type == "SBC_PML")` (~1489–1490) | **Delete branch**; keep `"SGBC"` path |
| `getSGBCTags()` | Will no longer pick up SBC_PML boundaries |

### `README.md`

| Section | Action |
|---------|--------|
| Boundary type `"SBC_PML"` | Remove |
| SBC_PML auto-generated layer docs (~169–171) | Remove |
| Cross-reference `docs/pml/` for volumetric PML | Add |

### Test data (review before delete)

| Path | Notes |
|------|-------|
| `testData/maxwellInputs/1D_SGBC_PML/` | Boundary SBC_PML — **remove or archive** |
| `testData/maxwellInputs/2D_Dipole_PML_G1/` | Likely boundary PML naming — **remove or repurpose** |
| `testData/maxwellInputs/2D_RCS_Circle_G1_PML/` | Likely boundary — **remove or repurpose** |
| `testData/maxwellInputs/2D_RCS_Circle_G2_PML/` | Likely boundary — **remove or repurpose** |

**Keep:** `testData/maxwellInputs/2D_RCS_Circle_Vol_PML/` (volumetric reference).

### `src/evolution/GlobalEvolution.cpp`

| Item | Action |
|------|--------|
| Comment "Boundary SGBC/SBC_PML" (~1013) | Update to SGBC only |

---

## Leftover volumetric PML experiments (user confirmed dead)

### JSON fields (not parsed today)

Present in `2D_RCS_Circle_Vol_PML.json` but **not implemented** in `assembleAttributeToMaterial()`:

- `"type": "PML"`
- `matches_vacuum`, `grading_order`, `target_reflection`, `kappa_max`, `alpha_max`, `active_axes`
- `"type": "vacuum"`

**Action:** Implement properly in A.2; not "remove" from JSON.

### `VolumetricRegionSubMesher`

| File | Status |
|------|--------|
| `src/components/SubMesher.h` | Class exists |
| `src/components/SubMesher.cpp` | Implementation + debug prints |

**Action for Milestone A:** **Keep code** but **do not wire** into solver unless needed for debug. Not required for physics. Optional: remove `[VolPML InterfaceDetect]` spam if unused.

**Not referenced** in `GlobalEvolution` or `Solver` today (grep confirms only SubMesher files).

---

## Keep unchanged (orthogonal features)

| Feature | Reason |
|---------|--------|
| `SGBC` boundary type | Still supported; independent of vol PML |
| `SGBCWrapper`, `SolverExtension` | Unchanged |
| TFSF / planewave / dipole | Unchanged |
| `buildSigmaPiecewiseVector()` | Bulk loss in non-PML materials |
| `collectGlobalConductiveOperator()` | Bulk loss; not PML |
| `external/mfem-geg/miniapps/dpg/util/pml.cpp` | MFEM fork utility; frequency-domain; **do not delete** |

---

## Hard-coded size audit (fix during A.4)

Grep targets when extending state:

```bash
rg "6 \* .*ndofs|6\*ndofs|6 \* fes" src/evolution src/solver src/components
```

Known locations:

- `GlobalEvolution.cpp`: ctor, Mult out size, ImplicitSolve assert
- `Fields.h`: SetSize
- `DGOperatorFactory.h`: globalRows (PML may stay 6N for curl block only; ψ separate)
- `ProbesManager.cpp`: ensure probes use 6N prefix only

---

## Git: `vol_pml` branch

User request:

1. Delete existing `vol_pml` branch (local and remote if present).
2. Branch fresh from **`dev`** (confirm exact base branch name at execution time).
3. Implement Milestone A on new `vol_pml`.

**Do not** reference old branch commits in implementation.

---

## README materials section (future)

After A.2, document:

```json
"type": "PML"
```

fields — see [03-json-schema-and-mesh.md](./03-json-schema-and-mesh.md).

Remove all `SBC_PML` boundary documentation.

---

## Verification after cleanup (A.1)

- [ ] `rg SBC_PML` returns no functional code (only docs/changelog if any)
- [ ] Project builds
- [ ] Existing SGBC case still runs (if present in testData)
- [ ] No parser references to `"type":"PML"` until A.2 (expected)
