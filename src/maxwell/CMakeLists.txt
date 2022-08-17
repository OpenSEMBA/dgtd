cmake_minimum_required(VERSION 3.0)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})


add_library(maxwell_solvers STATIC 
	"Solver.h" "Solver.cpp" 
	"FiniteElementEvolution.h" "FiniteElementEvolution.cpp" 
	"mfemExtension/BilinearIntegrators.h" "mfemExtension/BilinearIntegrators.cpp"
	"Types.h" 
	"Material.h" "Material.cpp"
	"Model.h" "Model.cpp" 
	"Probes.h" "Probes.cpp" 
	"Sources.h" "Sources.cpp"
)

target_link_libraries(maxwell_solvers mfem)

