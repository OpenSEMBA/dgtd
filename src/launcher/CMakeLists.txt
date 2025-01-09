message(STATUS "Creating build system for opensemba_dgtd")

add_executable(opensemba_dgtd "launcher.cpp")

target_link_libraries(opensemba_dgtd 
		maxwell-driver 
		mfem
)