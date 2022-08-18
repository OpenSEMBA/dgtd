cmake_minimum_required(VERSION 3.0)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(testHelpers 
	"GlobalFunctions.cpp" "GlobalFunctions.h"
 "AnalyticalFunctions1D.h" "AnalyticalFunctions2D.h")

target_link_libraries(testHelpers Eigen3::Eigen)

