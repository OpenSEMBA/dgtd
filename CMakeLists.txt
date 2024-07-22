cmake_minimum_required(VERSION 3.5)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd)

set(CMAKE_CXX_STANDARD 17)

option(MAXWELL_USE_MPI OFF)
option(MAXWELL_ENABLE_PYTHON_BINDINGS OFF)

include_directories(
	src/ 
	src/maxwell/
	src/adapter/
)

add_definitions(-DMFEM_USE_DOUBLE)
add_subdirectory(src)

if(MAXWELL_ENABLE_PYTHON_BINDINGS)
	add_subdirectory(pythonBindings)
endif()

enable_testing()

add_subdirectory(test/  )

add_test(maxwell ${CMAKE_BINARY_DIR}/bin/maxwell_tests)
add_test(mfem ${CMAKE_BINARY_DIR}/bin/mfem_tests)
add_test(hesthavenComparison ${CMAKE_BINARY_DIR}/bin/hesthavenComparison_tests)
add_test(adapter ${CMAKE_BINARY_DIR}/bin/adapter_tests)
add_test(rcs ${CMAKE_BINARY_DIR}/bin/rcs_tests)
