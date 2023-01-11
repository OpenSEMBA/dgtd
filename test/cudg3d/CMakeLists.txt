cmake_minimum_required(VERSION 3.0)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

include_directories(
	./
	${cudg3d_INCLUDE_DIRS}
)

message(STATUS "Creating build system for testMaxwell")
message(STATUS "Maxwell include dirs: " ${maxwell_INCLUDE_DIRS})
message(STATUS "MFEM include dirs: " ${MFEM_INCLUDE_DIRS})

add_definitions(-DPOLYNOMIAL_ORDER=3)

add_executable(cudg3d_tests 
	"dg/CellTest.cpp"
	"Cudg3dTest.cpp")



target_link_libraries(cudg3d_tests
	cudg3d
	opensemba_parsers
	mfem
	GTest::gtest GTest::gtest_main
)

