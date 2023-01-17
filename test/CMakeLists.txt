cmake_minimum_required(VERSION 3.0)

include_directories(
	./helpers
	${maxwell_INCLUDE_DIRS}
)

message(STATUS "Configuring build system for semba dgtd tests.")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_subdirectory(helpers)
add_subdirectory(maxwell)
add_subdirectory(mfem)
add_subdirectory(mfemHesthavenComparison)
add_subdirectory(thirdpartysoft)