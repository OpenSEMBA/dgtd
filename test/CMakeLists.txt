cmake_minimum_required(VERSION 3.0)

message(STATUS "Configuring build system for semba dgtd tests.")

add_subdirectory(jacobi)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(dgtd_tests ${SRCS})

target_link_libraries(dgtd_tests
                        GTest::gtest 
		                GTest::gtest_main
                        opensemba_core)

add_test(opensemba  ${CMAKE_BINARY_DIR}/bin/opensemba_tests)


