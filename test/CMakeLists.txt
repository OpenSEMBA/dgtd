cmake_minimum_required(VERSION 3.5)

include_directories(
	./
	${maxwell_INCLUDE_DIRS}
)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_subdirectory(maxwell)
add_subdirectory(mfem)
add_subdirectory(hesthavenComparison)
add_subdirectory(adapter)
add_subdirectory(rcs)