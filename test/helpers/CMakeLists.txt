cmake_minimum_required(VERSION 3.0)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})
include_directories(${maxwell_INCLUDE_DIRS})

message(STATUS "Creating build system for testHelpers")

add_library(testHelpers 
	"GlobalFunctions.cpp"
    )

target_link_libraries(testHelpers maxwell-solver Eigen3::Eigen GTest::gtest)

