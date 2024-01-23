cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for mfem_tests")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(FFTW3 CONFIG REQUIRED)

include_directories(
	./
)

add_executable(other_tests
"OtherTest.cpp"
)

target_link_libraries(other_tests 
	mfem
	GTest::gtest GTest::gtest_main
)

