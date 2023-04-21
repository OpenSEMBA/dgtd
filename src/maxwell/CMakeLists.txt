cmake_minimum_required(VERSION 3.0)

project(maxwell)

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_subdirectory(evolution)
add_subdirectory(mfemExtension)

add_library(maxwell STATIC 
	"Solver.cpp" 
	"Material.cpp"
	"Model.cpp" 
	"Sources.cpp"
    "ProbesManager.cpp" 
	"SourcesManager.cpp" 
	"Fields.cpp" 
 "PhysicalConstants.h")

target_link_libraries(maxwell maxwell-evolution mfemExtension mfem)
