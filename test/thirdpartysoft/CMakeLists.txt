cmake_minimum_required(VERSION 3.0)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

include_directories(
	./
	${maxwell_INCLUDE_DIRS}
)

add_executable(third_party_soft_tests
"GmshTest.cpp")

target_link_libraries(third_party_soft_tests mfem
	maxwell testHelpers
	GTest::gtest GTest::gtest_main
)


