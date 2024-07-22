cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for case_manager")

find_package(nlohmann_json CONFIG REQUIRED)

include_directories(./)

add_executable (cases_tests
	"CasesTest.cpp")

target_link_libraries(cases_tests maxwell-driver GTest::gtest GTest::gtest_main)