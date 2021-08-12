cmake_minimum_required(VERSION 3.0)

if(NOT DEFINED _VCPKG_INSTALLED_DIR AND DEFINED ENV{_VCPKG_INSTALLED_DIR})
   set(_VCPKG_INSTALLED_DIR $ENV{_VCPKG_INSTALLED_DIR})
endif()

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

macro(add_subdirectories dir bin)
    set(curdir ${CMAKE_CURRENT_LIST_DIR}/${dir})
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            if(EXISTS ${curdir}/${child}/CMakeLists.txt)
                add_subdirectory(${curdir}/${child} ${bin}/${dir}/${child})
            endif()
            if(EXISTS ${curdir}/${child}/src/CMakeLists.txt)
                add_subdirectory(${curdir}/${child}/src ${bin}/${dir}/${child})
            endif()
        endif()
    endforeach()
endmacro()

function(add_sources dir srcs)
    set(curdir ${CMAKE_CURRENT_LIST_DIR}/${dir})
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            if (dir STREQUAL ".")
                add_sources(${child} ${srcs})
            else()
                add_sources(${dir}/${child} ${srcs})
            endif()
        endif()
    endforeach()
    file(GLOB HDRS_LOCAL RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${curdir}/*.h)
    file(GLOB TMPL_LOCAL RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${curdir}/*.hpp)
    file(GLOB SRCS_LOCAL RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${curdir}/*.cpp)
    set(SRCS_LOCAL_ALL ${HDRS_LOCAL} ${TMPL_LOCAL} ${SRCS_LOCAL})
    set(${srcs} ${${srcs}} ${SRCS_LOCAL_ALL} PARENT_SCOPE)
    if (NOT dir STREQUAL ".")
        string(REPLACE "/" "\\\\" src_grp ${dir})
        source_group(${src_grp} FILES ${SRCS_LOCAL_ALL})
    endif()
endfunction()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(maxwell CXX)

if (CMAKE_COMPILER_IS_GNUCXX)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-std=c++11" compiler_supports_cxx11)
    check_cxx_compiler_flag("-std=c++0x" compiler_supports_cxx0X)
    if(compiler_supports_cxx11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    elseif(compiler_supports_cxx0X)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else()
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. "
                            "Please use a different C++ compiler.")
    endif()
endif()

add_subdirectory(src obj/)

include_directories(${CMAKE_CURRENT_LIST_DIR}/src)

add_custom_target(clobber WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                          COMMAND ${CMAKE_COMMAND} -E remove_directory obj
                          COMMAND ${CMAKE_COMMAND} -E remove_directory bin
                          COMMAND ${CMAKE_COMMAND} -E remove_directory lib)
			  

enable_testing()
add_subdirectory  (test/    obj/test/)
add_test(semba ${CMAKE_BINARY_DIR}/bin/maxwell_tests)
