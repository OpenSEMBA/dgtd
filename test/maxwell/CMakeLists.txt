cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for Maxwell_Tests")

message(STATUS "Maxwell include dirs: " ${maxwell_INCLUDE_DIRS})

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})
message(STATUS "MFEM include dirs: " ${MFEM_INCLUDE_DIRS})

add_executable(maxwell_tests
	"maxwell_tests.cpp"
	"Solver1DTest.cpp" 
	"Solver2DTest.cpp" 
	"Solver3DTest.cpp"  
	"MaterialTest.cpp" 
	"ProbesManagerTest.cpp" 
	"BilinearIntegratorsExtensionTest.cpp"
	"BilinearFormExtensionTest.cpp" 
	"LinearFormExtensionTest.cpp" "ParSolver2DTest.cpp")

target_link_libraries(maxwell_tests 
	maxwell mfem testHelpers 
	GTest::gtest
)

