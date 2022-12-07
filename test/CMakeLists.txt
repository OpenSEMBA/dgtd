cmake_minimum_required(VERSION 3.0)

message(STATUS "Configuring build system for semba dgtd tests.")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(dgtd_tests     
    "jacobi/RuleTest.cpp"
    "jacobi/LineTest.cpp"
    "jacobi/TriangleTest.cpp"
)

target_link_libraries(dgtd_tests dgtd_jacobi GTest::gtest GTest::gtest_main)

add_test(dgtd ${CMAKE_BINARY_DIR}/bin/dgtd_tests)


