cmake_minimum_required(VERSION 3.0)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

include_directories(${maxwell_solvers_INCLUDE_DIRS})

add_executable(maxwell_tests 
"TestMFEMFunctionality.cpp" 
"maxwell/TestSolver.cpp" 
"maxwell/TestBilinearIntegrators.cpp" 
"maxwell/TestMaterial.cpp" 
"maxwell/TestProbes.cpp" 
"maxwell/TestSources.cpp" 
"TestGlobalFunctions.cpp" "TestGlobalFunctions.h")

target_link_libraries(maxwell_tests maxwell_solvers GTest::gtest GTest::gtest_main)

add_executable(mfem_hesthaven_tests
"MFEMHesthavenComparison/TestMfemHesthaven1D.cpp" 
"MFEMHesthavenComparison/TestMfemHesthaven3D.cpp" 
"MFEMHesthavenComparison/TestMfemHesthavenFunctions.cpp" "MFEMHesthavenComparison/TestMfemHesthavenFunctions.h" 
"TestGlobalFunctions.cpp" "TestGlobalFunctions.h" 
"MFEMHesthavenComparison/TestMfemHesthavenTrace1D.cpp")

target_link_libraries(mfem_hesthaven_tests maxwell_solvers GTest::gtest GTest::gtest_main)


