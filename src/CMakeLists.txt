cmake_minimum_required(VERSION 3.0)

project(maxwell_solvers)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell_solvers STATIC 
	"Solver2D.cpp" "Solver2D.h"
	"maxwell1D/Solver.h" "maxwell1D/Solver.cpp" "maxwell1D/FE_Evolution.h" "maxwell1D/FE_Evolution.cpp" "maxwell1D/BilinearIntegrators.cpp" "maxwell1D/BilinearIntegrators.h"   )

