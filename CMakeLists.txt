cmake_minimum_required(VERSION 3.5)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd)

set(CMAKE_CXX_STANDARD 17)

option(MAXWELL_USE_MPI OFF)
option(MAXWELL_ENABLE_PYTHON_BINDINGS OFF)
option(MAXWELL_USE_MFEM_AS_SUBDIRECTORY ON)

include_directories(
	src/ 
	src/maxwell/
)

if (MAXWELL_USE_MFEM_AS_SUBDIRECTORY)
	add_subdirectory(external/mfem-geg)
endif()

add_definitions(-DMFEM_USE_DOUBLE)
add_subdirectory(src)

if(MAXWELL_ENABLE_PYTHON_BINDINGS)
	add_subdirectory(pythonBindings)
endif()

enable_testing()

add_subdirectory(test/  )

