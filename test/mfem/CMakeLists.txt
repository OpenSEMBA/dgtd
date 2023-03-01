cmake_minimum_required(VERSION 3.0)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

include_directories(
	./
)

add_executable(mfem_tests 
	"FiniteElementSpaceTest.cpp"
	"MeshTest.cpp" 
	"GridFunctionTest.cpp"
	"AlgebraTest.cpp")

target_link_libraries(mfem_tests 
	mfem 
	testHelpers
	GTest::gtest GTest::gtest_main
)

