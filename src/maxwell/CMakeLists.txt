cmake_minimum_required(VERSION 3.0)

project(maxwell)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell STATIC 
	"Solver.cpp" 
	"mfemExtension/BilinearIntegrators.cpp"
	"Material.cpp"
	"Model.cpp" 
	 
	"Sources.cpp"
    "ProbesManager.cpp" 
	"SourcesManager.cpp" 
	"Fields.cpp" 
	"MaxwellEvolution3D.cpp"
	"MaxwellEvolution2D.cpp"
	"MaxwellEvolution1D.cpp" 
	"MaxwellDefs.cpp"
	"MaxwellDefs1D.cpp"
	"mfemExtension/BilinearForm_IBFI.cpp" 
)

target_link_libraries(maxwell mfem Eigen3::Eigen)
