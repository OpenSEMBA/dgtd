# maxwell-dgtd
Maxwell's curl equations solver using discontinuous Galerkin methods


## Compiling

Compilation needs vcpkg with the following packages:
- eigen3
- gtest
- metis (for MPI)
- hypre (for MPI)

### OpenMP and MPI in Windows

- OpenMP requires an LLVM compiler. It has been tested Intel OneAPI (Base kit and HPC kit). To compile, use the following CMake Command Arguments when compiling MFEM and maxwell dgtd:
    ```
    -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icx
    ```
- MPI requires mfem to be compiled with MFEM_USE_MPI which requires Hypre and METIS.
  - Option 1, install with vcpkg. make sure to mark METIS_VERSION_5 and that th
  - Option 2, compiling from sources: Hypre can be cloned from https://github.com/hypre-space/hypre. It must be compiled **and installed** using a CMakeLists.txt in hypre/src. It is possible that the following variable has to be manually set:
    ```
        HYPRE_LIBRARIES=[PATH TO lib]\HYPRE.lib
    ```
    METIS 5.0.1 can be compiled in Windows with CMake. For VS2022+, comment out the line 
    ```
        #define rint(x) ((idx_t)((x)+0.5))
    ```
    in metis\GKlib\gk_arch.h

### OpenMP and MPI Linux

- Tested to compile with gcc.
- OpenMP has been tested to work.

### Defining a JSON file

Cases can be defined by parsing a JSON file with the problem information. See a complete [example](https://github.com/OpenSEMBA/dgtd/blob/alejandro/testData/maxwellInputs/1D_PEC_Centered/1D_PEC_Centered.json) of a valid JSON file. 
An OpenSemba/dgtd JSON file must have the following structure, **bold** entries are **required**:

- solver_options:
  	- Object. User can customise some of the simulation options in this section. If undefined, the problem will use the default settings.
	    - solver_type:
	        - String. Defines the evolution operator that will be used in the simulation. Can be "Centered" or "Upwind". If undefined, defaults to "Upwind".
	    - final_time:
	      	- Double. In Natural Units (1 meter/c), how long the problem will be simulated. If undefined, defaults to 2.0.
	    - time_step:
	      	- Double. In Natural Units (1 meter/c), defines the time step increments between iterations. Must be defined in 2D and/or 3D. Overrides CFL for 1D if both are defined. If undefined or set to 0.0 in 1D, will use an automatic time step calculation approach through CFL's condition.
	    - cfl:
	      	- Double. Courant–Friedrichs–Lewy condition, defines time step increments between iterations for 1D. Not available for 2D and/or 3D. If 'time_step' is defined this parameter is ignored. Defined between 0.0 and 1.0.
	  	- order:
	  	  	- Integer. Polynomial order of the interpolation function used in the Finite Element Collection. If undefined, defaults to 3.
	    - spectral:
	      	- Boolean. Use an Evolution Operator that uses a complete matrix form for all E and H unknowns but allows to calculate an approximate time step through analysis of the Eigenvalues of the matrix. Heavy computational cost, slower than the default Evolution Operators. Does not support the latest features. If undefined, defaults to false.
     
Example of a complete solver_options section:     
https://github.com/OpenSEMBA/dgtd/blob/5dd67ef6435066172e9387d674ea7ccc4c6a8b87/testData/maxwellInputs/1D_PEC_Centered/1D_PEC_Centered.json#L2-L8

- **Model**:
  	- Object. Contains geometrical information on the problem. User can define materials and boundaries in this section.
	  	- **filename**:
	  	  	- String. Name of the mesh file with the geometrical information for the problem. Must have the same name as the case folder.
	  	- **materials**:
	  	  	- Array of materials. Contains the physical information of the material domains defined in the problem. At least **one material must exist** in the problem and be defined. Each material is composed of the following parameters:
	  	  	  	- **tags**:
	  	  	  	  	- Array of integers. Contains the tags referred to the geometrical tags in the mesh that define the domain regions of the problem. (Volumes in 3D, Surfaces in 2D, Segments in 1D.)
	  	  	  	- **class**:
	  	  	  	  	- String. Defines the type of material.
	  	  	  	- relative_permittivity:
	  	  	  	  	- Double. Defines the relative permittivity parameter of the material. If undefined, defaults to 1.0.
	  	  	  	- relative_permeability:
	  	  	  	  	- Double. Defines the relative permeability parameter of the material. If undefined, defaults to 1.0.
	  	- **boundaries**:
	  	  	- Array of boundaries. Contains the physical information of the material interfaces or boundary faces defined in the problem. At least **one boundary must exist** in the problem and be defined. Each boundary is composed of the following parameters:
	  	  	  	- **tags**:
	  	  	  	  	- Array of integers. Contains the tags referred to the geometrical tags in the mesh that define the interface or boundary faces of the problem. (Surfaces in 3D, Segments in 2D, Points in 1D.)
	  	  	  	- **class**:
	  	  	  	  	- String. Defines the type of boundary to apply. Can be "PEC", "PMC", "SMA".
  	  	  	 
Example of a model section with two defined materials each with a different tag, Vacuum and a Dielectric with a specified permittivity; the model also has a boundary defined by two geometrical tags that share the same boundary type. The mesh and the folder for the case share the same name, with the exception of the mesh's file format:
https://github.com/OpenSEMBA/dgtd/blob/5dd67ef6435066172e9387d674ea7ccc4c6a8b87/testData/maxwellInputs/1D_PEC_Centered/1D_PEC_Centered.json#L11-L30

- Probes:
  	- Object. User can customise data extraction in this section. If undefined, no data extraction will be performed.
	  	- exporter:
	  	  	- Object. If defined, enables Paraview data exporting for posterior visualization.
	  	  	  	- steps:
	  	  	  	  	- Integer. Every how many time steps the solver will store data for visualization. All E and H componentts will be saved. If undefined, defaults to 1.
	  	- field:
	  	  	- Array. If defined, stores all fields at the specified point every time step.
	  	  	  	- position:
	  	  	  	  	- Array of doubles. Geometrical position at which data will be saved. Array must be n-Dimensional and define the required X, Y and/or Z coordinates. i.e. A 2D mesh only requires to define X and Y coordinates. ***Warning:*** If the point is defined outside the physical boundaries of the mesh, the simulation will crash.
	  	- surfaces:
	  	  	- Array. ***To be implemented.*** If defined, data will be extracted on the specified interfaces or boundary faces of the problem.
	  	  	  	- field:
	  	  	  	  	- String. Can be "E" or "H". Field to save.
	  	  	  	- tags:
	  	  	  	  	- Array of integers. Contains the tags referred to the geometrical tags in the mesh that define where to extract data.

Example of a probes section with all types of probes defined. An exporter probe that saves data at every step, a fields probe that stores all fields at position 0.0 for a 1D Mesh, and a surface probe that will store the "E" field at the surface with geometrical tag "1":
https://github.com/OpenSEMBA/dgtd/blob/2e153b3978c770c8fc5ef299de6bf6b51fd98ef0/testData/maxwellInputs/1D_PEC_Centered/1D_PEC_Centered.json#L32-L47

- **Sources**:
  	- Array of sources. User can customise problem illumination in this section. User should define at least one type of illumination from the available types. Different illumination types have different required parameters.
		- type:
  	  		- String. Can be "initial" or "totalField".
  	  	- tags:
  	  		- Array of integers. Geometrical tags for interfaces or boundary faces which define a totalField region. Only required if type: "totalField".
  	    - field_type:
  	      	- String. Only required if type: "initial". Can be "E" or "H".
  	    - center:
  	      	- Array of doubles. Center at with the illumination will be placed at the starting time. Only required if magnitude: "gaussian". Array must be n-Dimensional and define the required X, Y and/or Z coordinates. i.e. A 2D mesh only requires to define X and Y coordinates. If undefined, defaults to zero for each component.
  	    - polarization:
  	      	- Array of doubles. Direction of the field_type's polarization. Array must be 3D, all X, Y and Z vector magnitudes must be defined.
  	    - dimension:
  	      	- Integer. Which parameters to use in the (x+y+z)<sup>2</sup> term of the gaussian expression. Only required if magnitude: "gaussian".
  	    - propagation:
  	      	- Array of doubles. Array must be 3D, all X, Y and Z vector magnitudes must be defined. Only required if type: "totalField".
  	    - magnitude:
  	      	- Object. Defines paremeters relevant to the mathematical expression of the illumination. Not all parameters are required for specific magnitude types.
  	      	  	- type:
  	      	  	  	- String. Can be "gaussian" or "resonant". If type: "totalField", only "gaussian" is available.
  	      	  	- spread:
  	      	  	  	- Double. Defines the standard deviation \sigma of a [gaussian expression](https://wikimedia.org/api/rest_v1/media/math/render/svg/00cb9b2c9b866378626bcfa45c86a6de2f2b2e40). Only required if magnitude: "gaussian".
  	      	  	- delay:
  	      	  	  	- Double. Time delay for the gaussian to appear at the center coordinates (if center is not defined, at position (0.0, 0.0, 0.0)). Only required if magnitude: "gaussian". ***Warning:*** If the source type is "totalField", take into account the travelling wave will also have to physically travel to the boundaries to appear if the entry boundary is not located on at 0.0 for either component.
  	      	  	- modes:
  	      	  	  	- Array of integers. How many standing waves will fit in the specified geometry for each spatial component. Array must be n-Dimensional and define the required X, Y and/or Z modes. i.e. A 2D mesh only requires to define X and Y modes.

Example of a sources section with all available types of sources defined for documentation purposes. For most cases, a single source is enough:
- An initial gaussian source defined on E<sub>y</sub> at x = 0.5 (1D problem), the wave that will appear will be of type "gaussian" with the defined delay and spread.
- An initial resonant source on E<sub>y</sub> with 2 standing waves along the X axis (1D Problem).
- A totalField source defined on interfaces or boundary faces 1, 2, 3 and 4, polarised on the direction Y for the "E" field and with a propagation vector in the direction X, the wave that will appear will be of type "gaussian" with the defined delay and spread.

https://github.com/OpenSEMBA/dgtd/blob/5dd67ef6435066172e9387d674ea7ccc4c6a8b87/testData/maxwellInputs/1D_PEC_Centered/1D_PEC_Centered.json#L50-L83

  	
  	  
