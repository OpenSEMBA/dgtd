cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for maxwell-solver")

add_library(maxwell-solver STATIC 
	"Solver.cpp" 
    "ProbesManager.cpp" 
	"SourcesManager.cpp" 
)

target_link_libraries(maxwell-solver 
	maxwell-evolution 
	maxwell-components 
	maxwell-math
	mfem
)
