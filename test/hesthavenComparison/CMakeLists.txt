cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for hesthavenComparison_tests")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(hesthavenComparison_tests
	"Hesthaven1DTest.cpp" 
	"Hesthaven2DTest.cpp"
	"Hesthaven3DTest.cpp" 
	"HesthavenFunctions.cpp" 
)

target_link_libraries(hesthavenComparison_tests
	maxwell-solver
	Eigen3::Eigen
	GTest::gtest GTest::gtest_main
)


