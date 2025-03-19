message(STATUS "Creating build system for maxwell_Tests")

if (SEMBA_DGTD_ENABLE_MPI)
	set(MPI_SRCS "solver/ParSolver2DTest.cpp")
	add_compile_definitions(SEMBA_DGTD_ENABLE_MPI)
endif()

add_executable(maxwell_tests
	"driver/DriverTest.cpp"
	"components/MaterialTest.cpp" 
	"mfemExtension/BilinearIntegratorsExtensionTest.cpp"
	"mfemExtension/BilinearFormExtensionTest.cpp" 
	"mfemExtension/LinearFormExtensionTest.cpp" 
	"math/GeometryTest.cpp"
	"components/EigenvalueEstimator.cpp"
	"solver/ProbesManagerTest.cpp" 
	"solver/SourcesTest.cpp"
	"components/SubMesherTests.cpp" 
	"components/SphericalTest.cpp" 
	"math/FunctionTest.cpp")

add_executable(maxwell_solver_tests
	"MaxwellTestLauncher.cpp" 
	${MPI_SRCS}
	"solver/Solver1DTest.cpp" 
	"solver/Solver2DTest.cpp" 
	"solver/Solver3DTest.cpp")

if(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS)
	message(STATUS "Appending ExtensiveSolverTests to maxwell_solver_tests")
	target_sources(maxwell_solver_tests PRIVATE
		"solver/ExtensiveSolver2DTest.cpp")
endif()


target_link_libraries(maxwell_tests
	maxwell-solver
	maxwell-driver
	maxwell-components
	GTest::gtest
	GTest::gtest_main
)

target_link_libraries(maxwell_solver_tests
	maxwell-solver
	maxwell-driver
	GTest::gtest
)

include_directories("./")