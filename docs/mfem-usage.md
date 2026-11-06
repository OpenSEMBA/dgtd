# How this solver uses MFEM

The numerical results are the evidence that the time-domain Maxwell solver works. This note is about the layer underneath those results: what [MFEM](https://mfem.org) already provided, and where this repository took an MFEM object and specialized it. The build links the OpenSEMBA fork in [`external/mfem-geg`](../external/mfem-geg).

MFEM is a C++ library for high-order finite element methods on unstructured meshes, with MPI parallelism and a device programming model that covers CPU, OpenMP, and CUDA [1, 2]. Two lists follow.

**Used as provided.** The class already does the job. The solver passes an order, a mesh, a coefficient, or a device string and calls the method. `DG_FECollection` is the prototype: changing polynomial order is a constructor argument.

**Specialized.** MFEM has a related object, and this solver changes how it is used or what it contains. `SubMesh` is the prototype. The class can extract elements by attribute. `TotalFieldScatteredFieldSubMesher` decides which elements are total field, which are scattered field, and which faces carry the interface, then asks `SubMesh` to build the meshes.

A third kind of code is only named at the end. The Cartesian SC-PML auxiliary update, Debye, Lorentz, and the SGBC sub-step are the time-domain method. They call the objects below. They are not a modification of an MFEM class.

## Used as MFEM provides them

### Discontinuous basis functions

The solver constructs one scalar discontinuous collection for the whole run, in [`src/solver/Solver.cpp`](../src/solver/Solver.cpp):

```cpp
fec_{ opts_.evolution.order, model_.getMesh().Dimension(), opts_.basis_type},
fes_{ buildFiniteElementSpace(& model_.getMesh(), &fec_) },
```

`DG_FECollection` takes the polynomial order, the mesh dimension, and a basis family. The JSON field `basis_type` selects that family. The default is Gauss–Lobatto; Gauss–Legendre, Bernstein, and the uniform nodal families are the other documented choices ([`docs/json-input-format.md`](json-input-format.md)). Face and volume kernels ask the element for values or derivatives through `FiniteElement::CalcShape` and `CalcDShape`. They do not tabulate a basis.

A hand-written replacement is a nodal basis on the reference segment, triangle, quadrilateral, tetrahedron, and hexahedron, for every order the input file is allowed to request, in each of those families. The lasting part is the degree-of-freedom layout: which local index sits at which node, and how a face inherits a subset of those indices so that the mass matrix, the curl, the flux, the probes, and the halo all agree. Changing the order would then be a new implementation. The same collection is legal in one, two, and three dimensions. The CFL estimate builds a one-dimensional `DG_FECollection` of the same order only to read the reference-node spacing.

### The mesh bound to that space

`ParFiniteElementSpace` places the collection on a concrete mesh and numbers the degrees of freedom. `buildFiniteElementSpace` constructs it and exchanges face-neighbor tables, because a discontinuous method needs the degrees of freedom across every shared face ([`src/solver/Solver.cpp`](../src/solver/Solver.cpp)):

```cpp
auto fes = std::make_unique<mfem::ParFiniteElementSpace>(m, fec);
fes->ExchangeFaceNbrData();
```

`Mesh::LoadFromFile` ([`src/driver/driver.cpp`](../src/driver/driver.cpp)) supplies elements, faces, vertices, volumes, and attributes. `GetElementDofs` and `GetFE` are the maps from a cell to its basis nodes. Volume attributes become materials. Boundary attributes become boundary conditions. The fields object also builds a space with vector dimension 3 on the same collection ([`src/evolution/Fields.h`](../src/evolution/Fields.h)), so ParaView can see $(E_x, E_y, E_z)$ as one grid function. The operator itself is assembled on the scalar space.

Building that coupling by hand means a mesh database plus a degree-of-freedom map that is recomputed whenever the order or the basis changes, including ghost faces whose indices must match the neighbor. `ensureElementTypeIsSame` in [`src/components/Model.cpp`](../src/components/Model.cpp) rejects a parallel mesh that mixes element types. Tetrahedra and hexahedra are both supported, as separate meshes. MFEM can represent a mixed mesh; this solver's parallel path requires one element type.

### Geometry maps, normals, and quadrature tables

`ElementTransformation` maps a reference point to the physical cell and exposes the Jacobian and its weight. `FaceElementTransformations` is the same object for a face, with a fixed convention for element 1 and element 2, including a face owned by a neighboring rank. Quadrature points come from `IntRules.Get(geometry, order)`. `CalcOrtho` applied to the face Jacobian returns a vector along the normal whose length carries the surface measure. `Trans.Weight()` is the same geometric factor in scalar form.

On a straight simplex a hand-written cross product is short. On a high-order curved element the Jacobian varies inside the face, the measure is non-polynomial, and the sign is defined by MFEM's element-1 / element-2 orientation. An inverted normal flips every upwind flux. Reproducing that for tets and hexes, and for interior, boundary, and MPI faces, is a geometry library. `ElementTransformation::Transform` is the same map evaluated for a physical coordinate. Planewave sampling and PML profiles call it; those call sites are listed with the smaller stock uses below.

### The assembly loop

`ParBilinearForm` visits elements and faces, calls an integrator, and scatters a small dense matrix into a `SparseMatrix`. It understands attribute markers, interior faces, boundary faces, and shared faces whose columns land in the neighbor block. Domain terms use stock integrators unchanged: `DerivativeIntegrator` for a Cartesian derivative, `MassIntegrator` for a weighted mass matrix.

The hand-written alternative is that visit-and-scatter loop: owned elements, interior faces, boundary faces, ghost faces, per-attribute masks, and compressed-sparse-row insertion with the neighbor column offset. The Maxwell integrands that this loop calls are in the second part of the note. The loop itself is MFEM.

### Element mass inverse and piecewise coefficients

The strong-form scheme advances the fields with $M^{-1}$ times the spatial operator. On a discontinuous space the mass matrix is block diagonal, so the inverse is an element-local dense inversion. The factory builds it with stock classes ([`src/components/DGOperatorFactory.h`](../src/components/DGOperatorFactory.h)):

```cpp
Vector aux{pd_.model.buildEpsMuPiecewiseVector(f)};
PWConstCoefficient PWCoeff(aux);
auto res = std::make_unique<BF>(&fes_);
res->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));
```

`Model` fills a vector of $\varepsilon$ or $\mu$ indexed by attribute ([`src/components/Model.cpp`](../src/components/Model.cpp)). `PWConstCoefficient` and `InverseIntegrator` are MFEM. Conductivity uses the same mass path with a piecewise $\sigma$. The factory left-multiplies curl and flux blocks by this inverse through `mfem::Mult`.

A hand-written mass matrix repeats the quadrature of the basis products, then inverts one dense matrix per element. A global sparse factorization would discard the block-diagonal structure. Diagonal lumping changes the dispersion of a high-order Maxwell scheme. The element inverse keeps the consistent mass.

### Partitioned mesh and face-neighbor exchange

Each rank owns its elements. The numerical flux on a shared face needs the neighbor's trace. On rank 0 the driver calls `Mesh::GeneratePartitioning`, MFEM's METIS partition of the element graph ([`src/driver/driver.cpp`](../src/driver/driver.cpp)). `Model` builds the parallel mesh from that array ([`src/components/Model.cpp`](../src/components/Model.cpp)):

```cpp
serialMesh_ = Mesh(mesh);
ensureElementTypeIsSame(mesh);
pmesh_ = ParMesh(comm, serialMesh_, partitioning);
```

`ParFiniteElementSpace::num_face_nbr_dofs` is the width of the halo. The Maxwell operator is one local `SparseMatrix` with that halo stored as extra columns: $6N$ rows and $6(N+n_{\mathrm{nbr}})$ columns, where $N$ is the number of local scalar degrees of freedom per component. Every Runge–Kutta stage, `GlobalEvolution::Mult` calls `ParGridFunction::ExchangeFaceNbrData` ([`src/evolution/GlobalEvolution.cpp`](../src/evolution/GlobalEvolution.cpp)). MFEM packs the face degrees of freedom, exchanges them with the neighboring ranks, and writes `FaceNbrData()`.

Replacing this means a graph partitioner, a ghost layer, a pack format, and column indices that match the unpack order. MFEM can also build a `HypreParMatrix` for a global sparse solve and algebraic multigrid. This explicit scheme stores neighbor degrees of freedom as columns of a rectangular matrix. HYPRE is initialized in [`src/launcher/launcher.cpp`](../src/launcher/launcher.cpp) because the MPI build of MFEM expects it. The evolution algebra does not call it. The extra TFSF and SGBC constraint on that partition is in the second part.

### One device for CPU, OpenMP, and CUDA

The launcher selects `cpu`, `omp`, or `cuda` and constructs `mfem::Device` ([`src/launcher/launcher.cpp`](../src/launcher/launcher.cpp)). `mfem::Vector` residency goes through `UseDevice`, `Read`, `Write`, and `ReadWrite`. `mfem::forall` is the parallel loop for the active backend. `GlobalEvolution::GetMemoryClass` tells the ODE solver where stage vectors live. `Fields` stores $[E, H]$ and the auxiliaries in one device vector ([`src/evolution/Fields.h`](../src/evolution/Fields.h)).

A manual GPU port is a second copy of every pack and update kernel, plus explicit allocation and copy discipline at each MPI boundary. The same `Mult` source compiles for a CPU run and a CUDA run because of this memory manager. Two policies layered on top of it, the halo copy when MPI is not GPU-aware and the choice of sparse matvec, are in the second part.

### Time integrators

`GlobalEvolution` is an `mfem::TimeDependentOperator` with `Mult` and `ImplicitSolve`. `Solver::assignODESolver` selects a stock stepper from the JSON `ode_type`: `RK4Solver` by default, or backward Euler, the trapezoidal rule, implicit midpoint, or an SDIRK scheme ([`src/solver/Solver.cpp`](../src/solver/Solver.cpp)). The time loop is

```cpp
odeSolver_->Step(fields_.allDOFs(), time_, truedt);
```

`Step` allocates stage vectors in the memory class above and calls `Mult` or `ImplicitSolve` as the Butcher tableau requires. A hand-written RK4 is a small scalar routine and a memory bug on a device vector of fields plus auxiliaries. The implicit family adds stage coupling and a linear solve for the stage derivative $k$ in $y' = f(t, y)$.

Runs that include Cartesian PML, Debye, Lorentz, or SGBC must select RK4. The `Solver` constructor throws otherwise. Those extra terms are applied inside `Mult`. The nested SGBC problem is still a stock MFEM stepper: `buildSGBCSolverOptions` selects `SDIRK33Solver` ([`src/solver/SolverExtension.cpp`](../src/solver/SolverExtension.cpp)).

### Smaller stock uses

**ParaView.** `ParaViewDataCollection` writes the parallel VTK time series ([`src/solver/ProbesManager.cpp`](../src/solver/ProbesManager.cpp), [`src/launcher/mor2paraview.cpp`](../src/launcher/mor2paraview.cpp)), including high-order geometry subsampling. The manual replacement is a parallel VTK writer that understands discontinuous fields.

**Point probes.** `ParMesh::FindPoints` returns an element id and a reference `IntegrationPoint`. `GridFunction::GetValue` evaluates the field there ([`src/solver/ProbesManager.cpp`](../src/solver/ProbesManager.cpp)). The manual version is a spatial search and an inversion of the element map. There is no GSLIB call site under `src/`.

**Uniform refinement.** If the JSON model sets `refinement`, the driver calls `Mesh::UniformRefinement` before the finite-element space exists ([`src/driver/driver.cpp`](../src/driver/driver.cpp)). Attributes propagate to the children. The DG collection is constructed afterwards at the requested order.

**Curvature export.** `Mesh::SetCurvature` makes the VTK geometry match the solver order ([`src/solver/ProbesManager.cpp`](../src/solver/ProbesManager.cpp), [`src/components/RCSSurfaceExporter.cpp`](../src/components/RCSSurfaceExporter.cpp)). The call is stock. The quadrature-order policy for those curved elements is in the second part.

**Aliased views.** `Fields` allocates one `mfem::Vector` and `MakeRef`s each `ParGridFunction` onto a strided range ([`src/evolution/Fields.h`](../src/evolution/Fields.h)). Component $d$ of $E$ starts at index $dN$; component $d$ of $H$ starts at index $(d+3)N$. Probes and ParaView read the memory the stepper updates.

**Node coordinates.** `FiniteElement::GetNodes` plus `ElementTransformation::Transform` give the physical position of each degree of freedom ([`src/solver/SourcesManager.cpp`](../src/solver/SourcesManager.cpp)). PML profiles use `Transform` at quadrature points ([`src/components/PMLProfiles.cpp`](../src/components/PMLProfiles.cpp)). On a curved mesh this is the isoparametric map of whatever basis `basis_type` selected.

**Dense factorization.** `GlobalEvolution::ImplicitSolve` forms $J = I - \Delta t\, A$ and factors it with `DenseMatrixInverse` when the system is serial and small ([`src/evolution/GlobalEvolution.cpp`](../src/evolution/GlobalEvolution.cpp)). That branch is the one the SGBC sub-solver hits, using a two-slot cache so both the RK4 half-stage and full-stage sub-step sizes stay factored across the run. Larger systems on the same function use GMRES. The factorisation class is stock. It is separate from the element mass inverse above.

**Stock transfer.** `mfem::TransferMap` copies a grid function between a parent space and a submesh space. Probe extraction and the RCS surface exporter use it ([`src/solver/ProbesManager.cpp`](../src/solver/ProbesManager.cpp), [`src/components/RCSSurfaceExporter.cpp`](../src/components/RCSSurfaceExporter.cpp)). The signed add and subtract used by TFSF are in the second part.

## Specialized for this solver

These objects exist in MFEM, or exist as an extension point MFEM already documents (`Coefficient`, `BilinearFormIntegrator`, `SubMesh`, `IntRules`). The solver fills them in, changes their policy, or, in one case, adds a missing walk to the fork.

### Submeshes for TFSF, near-to-far, and PML regions

`SubMesh::CreateFromDomain` extracts the elements whose attributes are in a marker and keeps a parent-element id map. That is the whole stock operation. It does not know what a total-field region is.

`TotalFieldScatteredFieldSubMesher` ([`src/components/SubMesher.cpp`](../src/components/SubMesher.cpp)) copies the parent mesh and retags it before that call. For each interface face it reads `FaceElementTransformations`, decides which adjacent element is total field and which is scattered field, and writes `TotalFieldMarker` or `ScatteredFieldMarker` onto those elements. Dimension-specific passes do this in 1D, 2D, and 3D. `CreateFromDomain` then runs three times: a global interface neighborhood, a total-field mesh, and a scattered-field mesh. `createSubMeshFromParent` also copies boundary attributes onto the child and restores the original element attributes after extraction, because the temporary markers were only a selection key.

`NearToFarFieldSubMesher` does the same kind of retag on the elements that touch the integration surface, then extracts one `SubMesh`. `VolumetricRegionSubMesher` builds vacuum and PML markers from the JSON region tags, detects the vacuum–PML interface faces on the parent, and only then calls `CreateFromDomain`.

Around those meshes, [`src/components/SubMesher.h`](../src/components/SubMesher.h) computes orientations MFEM does not attach to a submesh: per-degree-of-freedom face normals (`computePerDofFaceNormals`), barycenter tests, and `buildFaceOrientation`. They call `CalcOrtho` and `Geometries.GetCenter`. The sign convention for the TFSF correction is this repository's.

A hand-written replacement is a region extractor, a parent-to-child element map, and those orientation tests, all consistent with the `DG_FECollection` numbering. `SubMesh` removes the extractor. It does not remove the decision of which elements belong to which region.

### Signed degree-of-freedom transfer

`MaxwellTransferMap` is a thin variant of `TransferMap`. The constructor still asks MFEM for the index table, through `SubMeshUtils::BuildVdofToVdofMap` and the submesh's `GetParentElementIDMap` ([`src/components/SubMesher.cpp`](../src/components/SubMesher.cpp)). `TransferAdd` and `TransferSub` then add or subtract along that table, which is what a TFSF update needs when the correction must accumulate. Stock `TransferMap` overwrites. The manual alternative is the same VDOF scatter, re-derived for every polynomial order.

### Maxwell face integrators

MFEM's assembly loop was listed above. The integrands are not stock DG trace integrators. `MaxwellDGZeroNormalJumpIntegrator`, `MaxwellDGOneNormalJumpIntegrator`, and `MaxwellDGTwoNormalJumpIntegrator` ([`src/mfemExtension/BilinearIntegrators.h`](../src/mfemExtension/BilinearIntegrators.h)) subclass `BilinearFormIntegrator` and implement `AssembleFaceMatrix`. They evaluate `CalcShape`, read the normal, and accumulate the penalty term, the single-normal coupling, and the product-of-normals coupling. `buildFaceMatrix` is the local outer product.

[`DGOperatorFactory`](../src/components/DGOperatorFactory.h) registers those integrators, calls `Assemble`, and places each scalar block into the $(E_x,\ldots,H_z)$ layout with the curl signs. `mergeBlocksToCSR` concatenates them into the rectangular operator. The visit over the mesh is MFEM. The component-wise flux and the block placement are this solver. Stock trace integrators do not encode that One / Zero / Two splitting, nor the PEC, PMC, and SMA coefficient tables.

### Near-to-far face integrals

`ParLinearForm` is the stock linear-form driver. `assembleLinearForm` in [`src/components/FarField.cpp`](../src/components/FarField.cpp) attaches `NTFFBdrFaceIntegrator` on the near-to-far marker. The integrator ([`src/mfemExtension/LinearIntegrators.cpp`](../src/mfemExtension/LinearIntegrators.cpp)) uses the same face transformations and `CalcOrtho` normals as the volume operator. The complex phase is a `FunctionCoefficient` lambda of the physical position. `FunctionCoefficient` itself is stock. The phase law and the integrator body are this repository's.

### A coefficient for the PML tensor

`PWConstCoefficient` and `FunctionCoefficient` are stock, and they were listed above where the solver only fills them. `SCPMLTensorCoefficient` is a subclass. Its `Eval` reads an `ElementTransformation` and an `IntegrationPoint` and returns one entry of the Bagci SC-PML tensor ([`src/components/DGOperatorFactory.h`](../src/components/DGOperatorFactory.h)). That object is passed to an ordinary `MassIntegrator` with an attribute marker. The assembly loop is unchanged. The tensor entries are the specialization. Writing a separate mass assembler for a spatially varying PML would repeat the quadrature MFEM already runs for every other coefficient.

### Interior boundaries in the fork

Interior PEC, PMC, SMA, and TFSF cuts are stored as boundary attributes, because the JSON tags them as boundaries. Both volume neighbors exist, so the flux needs a two-sided face. Stock assembly walks boundary faces as one-sided. The fork adds `Mesh::GetInternalBdrFaceTransformations` and `AddInternalBoundaryFaceIntegrator` on `BilinearForm` and `ParBilinearForm`, including shared faces in parallel. The factory calls that hook. The integrator still receives a `FaceElementTransformations`, two finite elements, and two degree-of-freedom arrays. The fork adds the walk over interior faces that carry a boundary attribute. It does not add a new mesh data structure. This is why the solver builds against `mfem-geg`.

### Quadrature order on curved elements

`IntRules` and `SetIntRule` are stock. The order this solver requests on curved meshes is not the library default. `buildDerivativeSubOperator` reads `GetNodalFESpace`. When the geometry order $q$ exceeds 1, it raises the `DerivativeIntegrator` rule by $(d-1)(q-1)$ on top of the degree-$p$ volume rule ([`src/components/DGOperatorFactory.h`](../src/components/DGOperatorFactory.h)). `setIntegrationRule` in [`src/mfemExtension/IntegratorFunctions.cpp`](../src/mfemExtension/IntegratorFunctions.cpp) adds a curvature term for faces. The reason is the summation-by-parts identity: a non-constant Jacobian increases the degree of the integrand, and a short rule shows up as a slow energy drift. The tables and the geometric factors remain MFEM's.

### Gmsh tags before the reader

`Mesh::LoadFromFile` reads a Gmsh 2.2 file and takes the physical tag as the MFEM attribute. `fixGmshMesh` runs first, on `.msh` files only ([`src/driver/driver.cpp`](../src/driver/driver.cpp)). In 3D it drops points and lines; in 2D it drops points. It copies the elementary tag into the physical-tag field (`tokens[3] = tokens[4]`) and renumbers the remaining elements. That adapted mesh is written to a temporary sidecar, which the reader loads and then deletes. The case `.msh` is left unchanged, so a thin-wire line is still there on the next launch. After that adaptation, the attribute-to-material maps agree with the JSON tags.

### Rank constraints after METIS

`GeneratePartitioning` balances the element graph. It has no notion that both sides of a TFSF or SGBC face must share a rank. `applyMetisPartitioningWithTFSFComponents` and the related pair-fix helpers in [`src/driver/driver.cpp`](../src/driver/driver.cpp) edit the partition array before `ParMesh` is built. The halo exchange then proceeds with the stock `ExchangeFaceNbrData` path. The primitive is MFEM's. The constraint is this solver's.

### Device policies on top of `mfem::Device`

`ExchangeFaceNbrData` and `SparseMatrix::Mult` are stock. Two call-site policies are not.

When MPI is not GPU-aware, the halo arrives in host memory. `sync_cuda_face_nbr_halos` in [`src/evolution/GlobalEvolution.cu`](../src/evolution/GlobalEvolution.cu) calls `Read()` so the device matvec sees it. The code branches on `Device::GetGPUAwareMPI()`.

cuSPARSE is turned off on every evolution matrix (`UseGPUSparse(false)` in `GlobalEvolution.cpp`). cuSPARSE records the output vector's size and requires it to equal the row count. The Maxwell block has $6N$ rows. The ODE state is taller once Cartesian PML auxiliaries are present ($12N$ in the code comment), and taller again with Debye and Lorentz. MFEM's CSR kernel, launched through `forall`, writes `height` rows and matches the host path. The matrix format stays MFEM's. The 2024 device paper also describes matrix-free partial assembly [2]. This solver stays on the assembled sparse matrix and applies it every stage.

## What is neither list

The two lists are MFEM objects and the specializations of those objects. The time-domain method sits on top of both.

`GlobalEvolution::Mult` applies the merged Maxwell matrix, then adds the TFSF source, the Cartesian SC-PML auxiliary update, Debye, Lorentz, and the SGBC face sub-step ([`src/evolution/GlobalEvolution.cpp`](../src/evolution/GlobalEvolution.cpp)). The ODE stepper sees the sum. Those updates are not a subclass of an MFEM integrator.

`HesthavenEvolution` is a deprecated reference operator. It uses an `L2_FECollection` on the same family of discontinuous spaces. `GlobalEvolution` is the path that carries PML, dispersive media, and SGBC.

## References

The papers below document the library. They are the source for the stock classes in the first list. They are not the source for the Maxwell fluxes, the TFSF submesh tagging, or the PML scheme.

1. R. Anderson, J. Andrej, A. Barker, J. Bramwell, J.-S. Camier, J. Cerveny, V. Dobrev, Y. Dudouit, A. Fisher, T. Kolev, W. Pazner, M. Stowell, V. Tomov, I. Akkerman, J. Dahm, D. Medina, and S. Zampini, "MFEM: A modular finite element methods library," *Computers & Mathematics with Applications*, vol. 81, pp. 42–74, 2021. doi: [10.1016/j.camwa.2020.06.009](https://doi.org/10.1016/j.camwa.2020.06.009).
2. J. Andrej, N. Atallah, J.-P. Bäcker, J. Camier, D. Copeland, V. Dobrev, Y. Dudouit, T. Duswald, B. Keith, D. Kim, T. Kolev, B. Lazarov, K. Mittal, W. Pazner, S. Petrides, S. Shiraiwa, M. Stowell, and V. Tomov, "High-performance finite elements with MFEM," *The International Journal of High Performance Computing Applications*, vol. 38, no. 5, pp. 447–467, 2024. doi: [10.1177/10943420241261981](https://doi.org/10.1177/10943420241261981).

The library itself is cited as MFEM, <https://mfem.org>, software identifier doi: [10.11578/DC.20171025.1248](https://doi.org/10.11578/DC.20171025.1248).
