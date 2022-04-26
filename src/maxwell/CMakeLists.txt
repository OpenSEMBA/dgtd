cmake_minimum_required(VERSION 3.0)

add_library(maxwell_solvers STATIC 
	"Solver1D.cpp" "Solver1D.h" 
	"FiniteElementEvolution.cpp" "FiniteElementEvolution.h" 
	"BilinearIntegrators.cpp" "BilinearIntegrators.h"    
	"Types.h" 
	 )

