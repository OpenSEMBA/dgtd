cmake_minimum_required(VERSION 3.25.2)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd)

set(CMAKE_CXX_STANDARD 17)

option(SEMBA_DGTD_ENABLE_MPI "Enable MPI support" OFF )
option(SEMBA_DGTD_ENABLE_OPEN_MP "Enable OpenMP support" ON )
option(SEMBA_DGTD_ENABLE_PYTHON_BINDINGS "Enable Python bindings" OFF )
option(SEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY "Use MFEM as a subdirectory" ON )

option(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS "Enable extensive case tests" ON )
option(SEMBA_DGTD_ENABLE_EXTENSIVE_SOLVER_TESTS "Enable extensive solver tests" ON )
option(SEMBA_DGTD_ENABLE_EXTENSIVE_RCS_TESTS "Enable extensive RCS tests" ON )

option(SEMBA_DGTD_ENABLE_TIMER_INFORMATION "Enable timer information" OFF)

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

add_subdirectory(src)

if(SEMBA_DGTD_ENABLE_PYTHON_BINDINGS)
	add_subdirectory(pythonBindings)
endif()

enable_testing()

add_subdirectory(test/)

