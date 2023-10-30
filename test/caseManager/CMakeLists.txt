cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for case_manager")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(nlohmann_json CONFIG REQUIRED)

add_library(case_manager
	"MaxwellCaseManager.cpp"
	"MaxwellModelManager.cpp"
	"MaxwellProbesManager.cpp"
	"MaxwellSourcesManager.cpp"
	"MaxwellSolverOptsManager.cpp")

target_link_libraries(case_manager
	mfem
	maxwell-solver
	nlohmann_json)