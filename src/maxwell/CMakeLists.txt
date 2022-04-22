cmake_minimum_required(VERSION 3.0)

add_library(maxwell_solvers STATIC 
	"Solver1D.cpp" "Solver1D.h" 
	"FE_Evolution.cpp" "FE_Evolution.h" 
	"BilinearIntegrators.cpp" "BilinearIntegrators.h"    
	"Types.h")

