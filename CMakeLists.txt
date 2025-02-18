cmake_minimum_required(VERSION 3.25.2)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd)

set(CMAKE_CXX_STANDARD 17)

option(MAXWELL_USE_MPI OFF)
option(MAXWELL_ENABLE_PYTHON_BINDINGS OFF)
option(MAXWELL_USE_MFEM_AS_SUBDIRECTORY ON)

option(MAXWELL_ENABLE_EXTENSIVE_CASE_TESTS ON)
option(MAXWELL_ENABLE_EXTENSIVE_SOLVER_TESTS ON)
option(MAXWELL_ENABLE_EXTENSIVE_RCS_TESTS ON)

option(MAXWELL_SHOW_TIMER_INFORMATION OFF)

include_directories(
	src/ 
	src/maxwell/
)

if (MAXWELL_USE_MFEM_AS_SUBDIRECTORY)
	add_subdirectory(external/mfem-geg)
else ()
	find_package(MFEM CONFIG REQUIRED)
	include_directories(${MFEM_INCLUDE_DIRS})
endif()

add_subdirectory(src)

if(MAXWELL_ENABLE_PYTHON_BINDINGS)
	add_subdirectory(pythonBindings)
endif()

if(MAXWELL_ENABLE_EXTENSIVE_CASE_TESTS)
	add_definitions(-DENABLE_EXTENSIVE_CASE_TESTS)
endif()

if(MAXWELL_ENABLE_EXTENSIVE_SOLVER_TESTS)
	add_definitions(-DENABLE_EXTENSIVE_SOLVER_TESTS)
endif()

if(MAXWELL_ENABLE_EXTENSIVE_RCS_TESTS)
	add_definitions(-DENABLE_EXTENSIVE_RCS_TESTS)
endif()

if(MAXWELL_SHOW_TIMER_INFORMATION)
	add_definitions(-DSHOW_TIMER_INFORMATION)
endif()

enable_testing()

add_subdirectory(test/  )

