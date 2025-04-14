cmake_minimum_required(VERSION 3.25.2)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd)

set(CMAKE_CXX_STANDARD 17)

option(SEMBA_DGTD_ENABLE_MPI "Enable MPI support" OFF )
option(SEMBA_DGTD_ENABLE_OPEN_MP "Enable OpenMP support" OFF )
option(SEMBA_DGTD_ENABLE_CUDA "Enable CUDA support" OFF )
option(SEMBA_DGTD_ENABLE_PYTHON_BINDINGS "Enable Python bindings" OFF )
option(SEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY "Use MFEM as a subdirectory" ON )

option(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS "Enable extensive case tests" ON )
option(SEMBA_DGTD_ENABLE_EXTENSIVE_SOLVER_TESTS "Enable extensive solver tests" ON )
option(SEMBA_DGTD_ENABLE_EXTENSIVE_RCS_TESTS "Enable extensive RCS tests" ON )

option(SEMBA_DGTD_ENABLE_TIMER_INFORMATION "Enable timer information" ON)

include_directories(
	src/ 
	src/maxwell/
)

if (SEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY)
	add_subdirectory(external/mfem-geg)
else ()
	find_package(MFEM CONFIG REQUIRED)
	include_directories(${MFEM_INCLUDE_DIRS})
endif()

if(SEMBA_DGTD_ENABLE_PYTHON_BINDINGS)
	add_subdirectory(pythonBindings)
endif()

macro(files_to_cuda EXE_SRCS)
    foreach(SRC_FILE IN LISTS ${EXE_SRCS})
        if (SEMBA_DGTD_ENABLE_CUDA)
            set_source_files_properties(${SRC_FILE} PROPERTIES LANGUAGE CUDA)
        endif()
    endforeach()
endmacro()

if(SEMBA_DGTD_ENABLE_CUDA)
    message(STATUS "Project files will be compiled with CUDA support")
    file(GLOB_RECURSE SRC_FILES "src/*.cpp" "test/*.cpp")
    files_to_cuda(SRC_FILES)
    
    if (NOT CMAKE_CUDA_HOST_COMPILER)
        set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    endif()
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD} CACHE STRING
        "CUDA standard to use.")
    set(CMAKE_CUDA_STANDARD_REQUIRED ON CACHE BOOL
        "Force the use of the chosen CUDA standard.")
    set(CMAKE_CUDA_EXTENSIONS OFF CACHE BOOL "Enable CUDA standard extensions.")
    set(CUDA_FLAGS "--expt-extended-lambda")
    if (CMAKE_VERSION VERSION_LESS 3.18.0)
        set(CUDA_FLAGS "-arch=${CUDA_ARCH} ${CUDA_FLAGS}")
    elseif (NOT CMAKE_CUDA_ARCHITECTURES)
        string(REGEX REPLACE "^sm_" "" ARCH_NUMBER "${CUDA_ARCH}")
        if ("${CUDA_ARCH}" STREQUAL "sm_${ARCH_NUMBER}")
        set(CMAKE_CUDA_ARCHITECTURES "${ARCH_NUMBER}")
        else()
        message(FATAL_ERROR "Unknown CUDA_ARCH: ${CUDA_ARCH}")
        endif()
    else()
        set(CUDA_ARCH "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")
    endif()
    message(STATUS "Using CUDA architecture: ${CUDA_ARCH}")
    if (CMAKE_VERSION VERSION_LESS 3.12.0)
        # CMake versions 3.8 and 3.9 require this to work; 3.10 and 3.11 are not
        # tested and may not actually need this (but should be ok to keep).
        set(CUDA_FLAGS "-ccbin=${CMAKE_CXX_COMPILER} ${CUDA_FLAGS}")
        set(CMAKE_CUDA_HOST_LINK_LAUNCHER ${CMAKE_CXX_COMPILER})
    endif()
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CUDA_FLAGS}")
    find_package(CUDAToolkit REQUIRED)
    set(CUSPARSE_FOUND TRUE)
    get_target_property(CUSPARSE_LIBRARIES CUDA::cusparse LOCATION)

endif()

add_subdirectory(src)

enable_testing()

add_subdirectory(test/)

