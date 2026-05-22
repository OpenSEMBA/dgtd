# Volumetric PML in dgtd

This document describes the first volumetric PML implementation in `dgtd`: what is implemented, how the main pieces interact, and what a user must provide to build a valid PML case.

## Scope

The current implementation is intentionally narrow:

- Only the `GlobalOperator` path is supported.
- The PML region must be an explicitly meshed outer shell.
- The shell must be vacuum-matched.
- The shell must be an axis-aligned wrapper around the physical region.
- The first formulation is sigma-only.
- `kappa = 1` and `alpha = 0` are kept internal.
- Serial and MPI execution are supported.

Not implemented in this first version:

- `HesthavenEvolution`
- `MaxwellEvolution`
- arbitrary PML shapes
- automatic shell generation
- user-exposed `kappa`, `alpha`, or `sigma_max`

## What Was Implemented

### Model and parsing

- `driver.cpp` recognizes `model.materials[]` entries with `"type": "PML"`.
- `Material` stays a simple scalar material holder. PML data is stored separately instead of creating a `PML_Material` subclass.
- `Model` stores PML metadata in `GeomTagToPMLRegion` and provides geometry helpers through `PMLBoxGeometry`.
- `Model::inferPMLBoxGeometry()` checks that the tagged region is a valid outer shell and computes inner and outer extents plus thickness per axis.

### Submesh and runtime state

- `VolumetricPMLSubMesher` builds the PML-only submesh using `SubMesh::CreateFromDomain` in serial and `ParSubMesh::CreateFromDomain` in MPI.
- `PMLWrapper` owns the PML submesh, the PML finite element space, the parent/submesh vdof mapping, the graded sigma profile, and the local correction operator.
- `PMLRegionState` stores the PML auxiliary state used across time steps.

### Operator path

- `DGOperatorFactory` can build a matched conductive operator on the PML submesh.
- The full-domain `buildGlobalOperator()` still provides the main Maxwell update.
- The PML operator is not a replacement for the global operator; it is a local correction applied only on the PML region.

### Time stepping

- `GlobalEvolution` builds and owns `PMLWrapper` when the model contains PML volumes.
- `Solver::step()` mirrors the SGBC pattern and now calls a PML checkpoint/finalize pair around the main ODE step.
- After the usual explicit update, `PMLWrapper` gathers the PML dofs from the parent fields, solves a local implicit correction on the PML submesh, and scatters the corrected values back into the parent fields.

### Diagnostics and protection

- Startup validation is handled in `driver.cpp`.
- The solver prints inferred thickness per active axis.
- When the source spectrum can be estimated, the solver also prints maximum frequency, shortest wavelength, and thickness-to-wavelength ratios.
- The solver prints normal PML cell size, cells across the thickness, effective dofs across the thickness, and the internal `sigma_max`.
- Invalid geometry hard-fails.
- Clearly underresolved shells hard-fail.
- Marginal shells warn.
- Very thick shells only report an informational message.

## How The Pieces Interact

The flow is:

1. `driver.cpp` parses a material entry tagged as `PML`.
2. `Model` stores the PML region metadata and infers the shell geometry.
3. `GlobalEvolution` creates a `PMLWrapper` when PML volumes are present.
4. `PMLWrapper` builds the PML submesh, submesh finite element space, sigma profile, and local matched conductive operator.
5. `Solver::step()` runs the usual explicit global step.
6. Right after that step, `PMLWrapper` applies an implicit correction only on the PML dofs.
7. The corrected PML values are written back into the parent field vectors.

In short: the main solver still advances the full Maxwell system as before, and the PML is added as a local post-step correction on a dedicated submesh.

## What The User Must Provide

To build a valid PML case, the user must provide all of the following:

- A mesh where the physical region and the PML shell are separate domain tags.
- A shell that wraps the physical region as an outer box.
- A shell that is axis-aligned.
- A shell tagged in `model.materials[]` with `"type": "PML"`.
- A physical region that remains a normal material entry, usually vacuum.

The user must not provide:

- manual PML thickness in JSON
- manual `sigma_max`
- manual `kappa` or `alpha`
- non-vacuum constitutive values on the PML material entry

The current JSON contract is based on `model.materials[]`, not on boundary conditions. A minimal pattern is:

```json
{
  "model": {
    "materials": [
      {
        "tags": [1],
        "type": "vacuum"
      },
      {
        "tags": [2],
        "type": "PML"
      }
    ]
  }
}
```

Here `1` is the physical region and `2` is the outer shell.

## Expectations For Mesh Design

When preparing a case, the user should expect these rules to matter:

- The shell must be the outermost region.
- The physical region must remain inside the shell with a clean box-like interface.
- The shell must exist in the mesh already; the solver will not create it.
- The shell should have enough cells across the normal thickness for the chosen polynomial order.
- If the source definition includes enough information to estimate frequency content, the solver will also judge thickness relative to wavelength.

Practical guidance:

- Use more than the bare minimum cell count across the shell.
- Keep the shell thickness consistent on opposite sides when possible.
- Start with vacuum physical material plus vacuum-matched PML shell.
- Prefer source descriptions with `spread` and, when relevant, `frequency`, so the wavelength-based diagnostics are available.

## Current Validation Surface

The implementation is covered by focused tests rather than a full production PML case mesh:

- parser and model tests for `type: PML`
- geometry inference tests
- serial and MPI submesh extraction tests
- direct `PMLWrapper` construction and damping tests
- a 1D end-to-end solver comparison showing stronger damping with the PML shell than with the same shell treated as ordinary vacuum

This means the core path is implemented and exercised, but a full external case still depends on the user supplying a purpose-built volumetric PML mesh.

## Summary

The current volumetric PML support is a real runtime feature, not only parser scaffolding. The solver can now:

- recognize PML regions from `model.materials[]`
- verify that the tagged region is a valid outer shell
- build a PML-only submesh
- assemble a local matched conductive correction
- apply an implicit PML update after the explicit global step
- protect the user with startup diagnostics before time stepping begins

To use it successfully, the user must provide a valid outer-shell mesh and tag that shell as a PML material region.