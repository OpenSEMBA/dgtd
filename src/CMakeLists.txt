if (SEMBA_DGTD_ENABLE_MPI)
	message(STATUS "SEMBA_DGTD_ENABLE_MPI:" ${SEMBA_DGTD_ENABLE_MPI})
	add_compile_definitions(SEMBA_DGTD_ENABLE_MPI)
endif()

add_subdirectory(launcher)
add_subdirectory(driver)
add_subdirectory(components)
add_subdirectory(evolution)
add_subdirectory(math)
add_subdirectory(mfemExtension)
add_subdirectory(solver)

