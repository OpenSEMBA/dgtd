cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for hesthavenComparison_tests")

add_executable(hesthavenComparison_tests
	"HesthavenFunctions.cpp" 
	"Hesthaven1DTest.cpp" 
	"Hesthaven2DTest.cpp"
	"Hesthaven3DTest.cpp" 
)

target_link_libraries(hesthavenComparison_tests
	maxwell-solver
	testHelpers
	GTest::gtest GTest::gtest_main
)


