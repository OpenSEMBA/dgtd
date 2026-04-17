message(STATUS "Creating build system for maxwell_Tests")

if (SEMBA_DGTD_ENABLE_MPI)
	add_compile_definitions(SEMBA_DGTD_ENABLE_MPI)
endif()

add_executable(maxwell_tests
	$<TARGET_OBJECTS:TestLauncher>
	"driver/DriverTest.cpp"
	"components/EigenvalueEstimator.cpp"
	"components/MaterialTest.cpp" 
	"components/ModelTest.cpp"
	"components/SIBCPoleResidueFittingTest.cpp"
	"components/SphericalTest.cpp" 
	"components/SubMesherTests.cpp" 
	"components/DGOperatorHelpersTest.cpp"
	"components/FarFieldTest.cpp"
	"mfemExtension/BilinearIntegratorsExtensionTest.cpp"
	"mfemExtension/BilinearFormExtensionTest.cpp" 
	"mfemExtension/LinearFormExtensionTest.cpp"
	"mfemExtension/IntegratorFunctionsTest.cpp"
	"math/GeometryTest.cpp"
	"math/PhysicalConstantsTest.cpp"
	"math/CalculusTest.cpp"
	"math/EigenMfemToolsTest.cpp"
	"solver/ProbesManagerTest.cpp" 
	"solver/SourcesTest.cpp"
	"math/FunctionTest.cpp")

add_executable(maxwell_solver_tests
	$<TARGET_OBJECTS:TestLauncher>
	${MPI_SRCS}
	"solver/Solver1DTest.cpp" 
	"solver/Solver2DTest.cpp" 
	"solver/Solver3DTest.cpp")

if(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS)
	message(STATUS "Appending Extensive Tests to maxwell_solver_tests")
	target_sources(maxwell_solver_tests PRIVATE
		"solver/ExtensiveSolver2DTest.cpp"
		"solver/ExtensionTest.cpp")
endif()


target_link_libraries(maxwell_tests
	maxwell-solver
	maxwell-driver
	maxwell-components
	GTest::gtest
)

target_link_libraries(maxwell_solver_tests
	maxwell-solver
	maxwell-driver
	GTest::gtest
)

include_directories("./")