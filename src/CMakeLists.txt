cmake_minimum_required(VERSION 3.0)

if (MAXWELL_USE_MPI)
	message(STATUS "MAXWELL_USE_MPI:" ${MAXWELL_USE_MPI})
	add_compile_definitions(MAXWELL_USE_MPI)
endif()

#add_subdirectory(cudg3d)
add_subdirectory(maxwell)
