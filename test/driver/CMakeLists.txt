cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for case_manager")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(nlohmann_json CONFIG REQUIRED)

include_directories(./)

add_executable (adapter_tests
	"AdapterTest.cpp")

target_link_libraries(adapter_tests
	adapter
	GTest::gtest GTest::gtest_main)