cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for maxwell_Tests")

if (MAXWELL_USE_MPI)
	set(MPI_SRCS "solver/ParSolver2DTest.cpp")
	add_compile_definitions(MAXWELL_USE_MPI)
	message(STATUS "MPI_SRCS=" ${MPI_SRCS})
endif()

add_executable(maxwell_tests
	"MaxwellSystemTestLauncher.cpp"
	${MPI_SRCS}
	"solver/Solver1DTest.cpp" 
	"solver/Solver2DTest.cpp" 
	"solver/Solver3DTest.cpp"  
	"components/MaterialTest.cpp" 
	"solver/ProbesManagerTest.cpp" 
	"solver/SourcesTest.cpp"
	"mfemExtension/BilinearIntegratorsExtensionTest.cpp"
	"mfemExtension/BilinearFormExtensionTest.cpp" 
	"mfemExtension/LinearFormExtensionTest.cpp" 
)

target_link_libraries(maxwell_tests 
	maxwell-solver
	testHelpers 
	GTest::gtest
)

