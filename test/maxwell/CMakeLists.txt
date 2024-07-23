cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for maxwell_Tests")

if (MAXWELL_USE_MPI)
	set(MPI_SRCS "solver/ParSolver2DTest.cpp")
	add_compile_definitions(MAXWELL_USE_MPI)
endif()

add_executable(maxwell_tests
	"driver/DriverTest.cpp"
	"components/MaterialTest.cpp" 
	"mfemExtension/BilinearIntegratorsExtensionTest.cpp"
	"mfemExtension/BilinearFormExtensionTest.cpp" 
	"mfemExtension/LinearFormExtensionTest.cpp" 
	"math/GeometryTest.cpp"
	"components/EigenvalueEstimator.cpp"
)

add_executable(maxwell_solver_tests
	"MaxwellTestLauncher.cpp" 
	${MPI_SRCS}
	"solver/Solver1DTest.cpp" 
	"solver/Solver2DTest.cpp" 
	"solver/Solver3DTest.cpp"  
	"solver/ProbesManagerTest.cpp" 
	"solver/SourcesTest.cpp"
)

target_link_libraries(maxwell_tests
	maxwell-solver
	maxwell-driver
	GTest::gtest
	GTest::gtest_main
)

target_link_libraries(maxwell_solver_tests
	maxwell-solver
	maxwell-driver
	GTest::gtest
)

include_directories("./")