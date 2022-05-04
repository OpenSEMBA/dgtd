cmake_minimum_required(VERSION 3.0)

add_library(maxwell_solvers STATIC 
	"Solver.cpp" "Solver.h" 
	"FiniteElementEvolution.cpp" "FiniteElementEvolution.h" 
	"BilinearIntegrators.cpp" "BilinearIntegrators.h"    
	"Types.h" "Options.h"
	"Material.h" "Material.cpp"
	"Model.h" "Model.cpp" 
	"Probes.h" "Probes.cpp" 
	"Sources.h" "Sources.cpp"  )

