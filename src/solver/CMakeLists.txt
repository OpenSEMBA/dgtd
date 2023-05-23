cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for maxwell-solver")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell-solver STATIC 
	"Solver.cpp" 
    "ProbesManager.cpp" 
	"SourcesManager.cpp" 
)

target_link_libraries(maxwell-solver maxwell-evolution maxwell-components mfem)
