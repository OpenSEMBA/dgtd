cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for maxwell_Tests")

if (MAXWELL_USE_MPI)
	set(MPI_SRCS "ParSolver2DTest.cpp")
	add_compile_definitions(MAXWELL_USE_MPI)
	message(STATUS "MPI_SRCS=" ${MPI_SRCS})
endif()

add_executable(maxwell_tests
	"MaxwellSystemTestLauncher.cpp"
	${MPI_SRCS}
	"Solver1DTest.cpp" 
	"Solver2DTest.cpp" 
	"Solver3DTest.cpp"  
	"MaterialTest.cpp" 
	"ProbesManagerTest.cpp" 
	"SourcesTest.cpp"
	"BilinearIntegratorsExtensionTest.cpp"
	"BilinearFormExtensionTest.cpp" 
	"LinearFormExtensionTest.cpp" 
)

target_link_libraries(maxwell_tests 
	maxwell-solver
	testHelpers 
	GTest::gtest
)

