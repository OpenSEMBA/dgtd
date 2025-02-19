cmake_minimum_required(VERSION 3.25.2)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd)

set(CMAKE_CXX_STANDARD 17)

option(SEMBA_DGTD_ENABLE_MPI OFF "Enable MPI support")
option(SEMBA_DGTD_ENABLE_PYTHON_BINDINGS OFF "Enable Python bindings")
option(SEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY ON "Use MFEM as a subdirectory")

option(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS ON "Enable extensive case tests")
option(SEMBA_DGTD_ENABLE_EXTENSIVE_SOLVER_TESTS ON "Enable extensive solver tests")
option(SEMBA_DGTD_ENABLE_EXTENSIVE_RCS_TESTS ON "Enable extensive RCS tests")

option(SEMBA_DGTD_ENABLE_TIMER_INFORMATION OFF "Enable timer information")

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

if(SEMBA_DGTD_ENABLE_EXTENSIVE_SOLVER_TESTS)
	add_definitions(-DENABLE_EXTENSIVE_SOLVER_TESTS)
endif()

if(SEMBA_DGTD_ENABLE_EXTENSIVE_RCS_TESTS)
	add_definitions(-DENABLE_EXTENSIVE_RCS_TESTS)
endif()


enable_testing()

add_subdirectory(test/  )

