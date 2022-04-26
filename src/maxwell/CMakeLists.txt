cmake_minimum_required(VERSION 3.0)

add_library(maxwell_solvers STATIC 
	"Solver1D.cpp" "Solver1D.h" 
	"FiniteElementEvolutionSimple.cpp" "FiniteElementEvolutionSimple.h" 
	"BilinearIntegrators.cpp" "BilinearIntegrators.h"    
	"Types.h" 
	"FiniteElementEvolutionMaterials.cpp" "FiniteElementEvolutionMaterials.h")

