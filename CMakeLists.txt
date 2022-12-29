cmake_minimum_required(VERSION 3.0)

if(NOT DEFINED _VCPKG_INSTALLED_DIR AND DEFINED ENV{_VCPKG_INSTALLED_DIR})
   set(_VCPKG_INSTALLED_DIR $ENV{_VCPKG_INSTALLED_DIR})
endif()
set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd CXX)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_subdirectory(external/core/)
include_directories(src/ external/core/src/ external/core/src/core/)

add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-DPOLYNOMIAL_ORDER=3)

add_subdirectory(src)

include_directories(src/cudg3d/)	  

enable_testing()

add_subdirectory(test/  )

add_test(maxwell ${CMAKE_BINARY_DIR}/bin/maxwell_tests)
add_test(mfem ${CMAKE_BINARY_DIR}/bin/mfem_tests)
add_test(mfem_hesthaven ${CMAKE_BINARY_DIR}/bin/mfem_hesthaven_tests)
add_test(3PTests ${CMAKE_BINARY_DIR}/bin/third_party_soft_tests)
