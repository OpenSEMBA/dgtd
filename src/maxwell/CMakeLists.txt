cmake_minimum_required(VERSION 3.0)

project(maxwell)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell STATIC 
	"SolverOptions.h" 		
	"Solver.h" "Solver.cpp" 
	"FiniteElementEvolution.h" "FiniteElementEvolution.cpp" 
	"mfemExtension/BilinearIntegrators.h" "mfemExtension/BilinearIntegrators.cpp"
	"Types.h" 
	"Material.h" "Material.cpp"
	"Model.h" "Model.cpp" 
	"Probes.h" "Probes.cpp" 
	"Sources.h" "Sources.cpp"
    "ProbesManager.h" "ProbesManager.cpp" 
	"SourcesManager.h" "SourcesManager.cpp" 
	"Fields.h" "Fields.cpp" 
)

target_link_libraries(maxwell mfem)
