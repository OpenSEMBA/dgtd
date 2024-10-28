cmake_minimum_required(VERSION 3.5)

if (MAXWELL_USE_MPI)
	message(STATUS "MAXWELL_USE_MPI:" ${MAXWELL_USE_MPI})
	add_compile_definitions(MAXWELL_USE_MPI)
endif()

add_subdirectory(launcher)
add_subdirectory(driver)
add_subdirectory(components)
add_subdirectory(evolution)
add_subdirectory(math)
add_subdirectory(mfemExtension)
add_subdirectory(solver)

