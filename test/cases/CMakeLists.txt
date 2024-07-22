cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for case_manager")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(nlohmann_json CONFIG REQUIRED)

include_directories(./)

add_executable (cases_tests
	"casesTest.cpp")

target_link_libraries(cases_tests maxwell-driver GTest::gtest GTest::gtest_main)