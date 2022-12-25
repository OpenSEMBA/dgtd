cmake_minimum_required(VERSION 3.0)

message(STATUS "Configuring build system for semba dgtd tests.")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(dgtd_tests
 "Cudg3dTest.cpp"
 "MFEMTest.cpp")

target_link_libraries(dgtd_tests 
	dgtd_solver
	dgtd_dg
	opensemba_parsers
	mfem
	GTest::gtest GTest::gtest_main
)

add_test(dgtd ${CMAKE_BINARY_DIR}/bin/dgtd_tests)


