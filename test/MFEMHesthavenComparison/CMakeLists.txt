cmake_minimum_required(VERSION 3.0)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(mfem_hesthaven_tests
	"MfemHesthavenFunctionsTest.cpp" 
	"MfemHesthaven1DTest.cpp" 
#	"MfemHesthaven3DTest.cpp" 
#	"MfemHesthavenTrace1DTest.cpp" 
#	"MfemHesthaven2DTest.cpp"
)

target_link_libraries(mfem_hesthaven_tests mfem
	maxwell testHelpers
	GTest::gtest GTest::gtest_main
)


