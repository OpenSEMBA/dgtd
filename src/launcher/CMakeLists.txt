message(STATUS "Creating build system for opensemba_dgtd")

add_executable(opensemba_dgtd "launcher.cpp")

target_link_libraries(opensemba_dgtd 
		maxwell-driver 
		mfem
)

if (SEMBA_DGTD_ENABLE_OPEN_MP)
	message(STATUS "Linking opensemba_dgtd with OpenMP libraries")
	target_link_libraries(opensemba_dgtd OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()