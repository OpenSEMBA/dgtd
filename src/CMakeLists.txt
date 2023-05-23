cmake_minimum_required(VERSION 3.0)

project(maxwell)

if (MAXWELL_USE_MPI)
	message(STATUS "MAXWELL_USE_MPI:" ${MAXWELL_USE_MPI})
	add_compile_definitions(MAXWELL_USE_MPI)
endif()


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
	"ProblemDescription.cpp"
)

target_link_libraries(maxwell maxwell-evolution mfemExtension mfem)
