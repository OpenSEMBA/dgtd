cmake_minimum_required(VERSION 3.0)

project(maxwell)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell STATIC 
	"Solver.cpp" 
	"Material.cpp"
	"Model.cpp" 
	"Sources.cpp"
    "ProbesManager.cpp" 
	"SourcesManager.cpp" 
	"Fields.cpp" 
	"Evolution3D.cpp"
	"Evolution2D.cpp"
	"Evolution1D.cpp" 
	"Evolution1D_Spectral.cpp"
	"Evolution2D_Spectral.cpp" 
	"Evolution3D_Spectral.cpp"
	"EvolutionMethods.cpp"
	"mfemExtension/BilinearForm_IBFI.cpp" 
	"mfemExtension/LinearForm_IBFI.cpp" 
	"mfemExtension/BilinearIntegrators.cpp"
	"mfemExtension/LinearIntegrators.cpp" 
	"mfemExtension/IntegratorFunctions.cpp" 
	"Calculus.h"
)

target_link_libraries(maxwell mfem Eigen3::Eigen)
