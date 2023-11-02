cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for case_manager")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(nlohmann_json CONFIG REQUIRED)

add_library(adapter
	"MaxwellAdapter.cpp"
	"ModelAdapter.cpp"
	"ProbesAdapter.cpp"
	"SourcesAdapter.cpp"
	"SolverOptsAdapter.cpp")

add_executable (adapter_tests
	"AdapterTest.cpp")

target_link_libraries(adapter
	mfem
	maxwell-solver
	nlohmann_json)

target_link_libraries(adapter_tests
	adapter
	mfem
	GTest::gtest GTest::gtest_main
	maxwell-solver
	nlohmann_json)