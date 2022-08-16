cmake_minimum_required(VERSION 3.0)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

include_directories(
	./
	${maxwell_solvers_INCLUDE_DIRS}
)

add_executable(maxwell_tests 
	"TestSolver.cpp" 
	"TestBilinearIntegrators.cpp" 
	"TestMaterial.cpp" 
	"TestProbes.cpp" 
	"TestSources.cpp"
	"TestMFEMFunctionality.cpp"
)

target_link_libraries(maxwell_tests 
	maxwell_solvers testHelpers 
	GTest::gtest GTest::gtest_main
)

