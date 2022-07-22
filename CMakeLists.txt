cmake_minimum_required(VERSION 3.0)

if(NOT DEFINED _VCPKG_INSTALLED_DIR AND DEFINED ENV{_VCPKG_INSTALLED_DIR})
   set(_VCPKG_INSTALLED_DIR $ENV{_VCPKG_INSTALLED_DIR})
endif()

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(maxwell)

add_subdirectory(src)

include_directories(${CMAKE_CURRENT_LIST_DIR}/src)	  

enable_testing()
add_subdirectory  (test/  )
add_test(maxwell ${CMAKE_BINARY_DIR}/bin/maxwell_tests)
add_test(mfem_hesthaven ${CMAKE_BINARY_DIR}/bin/mfem_hesthaven_tests)
