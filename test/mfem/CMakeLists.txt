cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for mfem_tests")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(
	./
)

add_executable(mfem_tests 
	"FiniteElementSpaceTest.cpp"
	"MeshTest.cpp" 
	"GridFunctionTest.cpp"
	"AlgebraTest.cpp"
)

target_link_libraries(mfem_tests 
	mfem 
	testHelpers
	Eigen3::Eigen
	GTest::gtest GTest::gtest_main
)

