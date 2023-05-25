cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for maxwell_Tests")

if (MAXWELL_USE_MPI)
	set(MPI_SRCS "solver/ParSolver2DTest.cpp")
	add_compile_definitions(MAXWELL_USE_MPI)
	message(STATUS "MPI_SRCS=" ${MPI_SRCS})
endif()

option(MAXWELL_TEST_ALLOW_PARAVIEW_EXPORT OFF)
if (MAXWELL_TEST_ALLOW_PARAVIEW_EXPORT)
	add_compile_definitions(MAXWELL_TEST_ALLOW_PARAVIEW_EXPORT)
endif()

add_executable(maxwell_tests
	"MaxwellTestLauncher.cpp"
	"components/MaterialTest.cpp" 
	"mfemExtension/BilinearIntegratorsExtensionTest.cpp"
	"mfemExtension/BilinearFormExtensionTest.cpp" 
	"mfemExtension/LinearFormExtensionTest.cpp" 
	${MPI_SRCS}
	"solver/Solver1DTest.cpp" 
	"solver/Solver2DTest.cpp" 
	"solver/Solver3DTest.cpp"  
	"solver/ProbesManagerTest.cpp" 
	"solver/SourcesTest.cpp"
	"adapter/OpensembaAdapterTest.cpp"
)

include_directories(./)

target_link_libraries(maxwell_tests 
	maxwell-solver
	GTest::gtest
)

