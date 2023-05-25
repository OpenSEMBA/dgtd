cmake_minimum_required(VERSION 3.0)

include_directories(
	./helpers
	${maxwell_INCLUDE_DIRS}
)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_subdirectory(helpers)
add_subdirectory(maxwell)
add_subdirectory(mfem)
add_subdirectory(hesthavenComparison)