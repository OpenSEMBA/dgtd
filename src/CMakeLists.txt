cmake_minimum_required(VERSION 3.5)

if (MAXWELL_USE_MPI)
	message(STATUS "MAXWELL_USE_MPI:" ${MAXWELL_USE_MPI})
	add_compile_definitions(MAXWELL_USE_MPI)
endif()

if (NOT MAXWELL_USE_MFEM_AS_SUBDIRECTORY)
	find_package(MFEM CONFIG REQUIRED)
	include_directories(${MFEM_INCLUDE_DIRS})
endif()

add_subdirectory(driver)
add_subdirectory(components)
add_subdirectory(evolution)
add_subdirectory(math)
add_subdirectory(mfemExtension)
add_subdirectory(solver)
