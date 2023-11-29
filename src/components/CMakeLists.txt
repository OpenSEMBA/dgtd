cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for maxwell-components")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell-components STATIC 
	"Material.cpp"
	"Model.cpp" 
	"Sources.cpp"
	"../solver/SolverInput.h"
	"EigenvalueEstimator.cpp"
	"OptimizationManager.cpp"
	"SubMesher.cpp"
	)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(maxwell-components mfem Eigen3::Eigen)
