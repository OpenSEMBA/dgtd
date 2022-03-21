cmake_minimum_required(VERSION 3.0)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

include_directories(${maxwell_solvers_INCLUDE_DIRS})

add_executable(maxwell_tests "TestDG.cpp" "TestSolver2D.cpp" "maxwell1D/TestSolver.cpp" "maxwell1D/TestBilinearIntegrators.cpp")

target_link_libraries(maxwell_tests maxwell_solvers mfem GTest::gtest GTest::gtest_main) 


