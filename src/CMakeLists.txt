cmake_minimum_required(VERSION 3.0)

project(maxwell_solvers)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell_solvers STATIC 
	"Solver.cpp" "Solver.h"
	"Solver1D.cpp" "Solver1D.h"
	"Solver1DRunge.cpp" "Solver1DRunge.h")

