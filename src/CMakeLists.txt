cmake_minimum_required(VERSION 3.0)



if (MAXWELL_USE_MPI)
	message(STATUS "MAXWELL_USE_MPI:" ${MAXWELL_USE_MPI})
	add_compile_definitions(MAXWELL_USE_MPI)
endif()

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_subdirectory(components)
add_subdirectory(adapter)
add_subdirectory(mfemExtension)
add_subdirectory(evolution)
add_subdirectory(solver)
