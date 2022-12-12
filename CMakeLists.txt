cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
    set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-DPOLYNOMIAL_ORDER=3)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd CXX)
add_subdirectory(external/core/)
include_directories(src/ external/core/src/ external/core/src/core/)

option(SEMBA_DGTD_USE_MPI OFF)

message(STATUS "Creating build system for dgtd with the following options:")
message(STATUS "SEMBA_DGTD_USE_MPI: " ${SEMBA_DGTD_USE_MPI})

add_subdirectory(src/)

enable_testing()
add_subdirectory(test/)