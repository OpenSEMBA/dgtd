message(STATUS "Creating build system for opensemba_dgtd")

find_package(PkgConfig REQUIRED)
pkg_check_modules(GSL REQUIRED gsl)

add_executable(opensemba_dgtd "launcher.cpp")

target_include_directories(opensemba_dgtd PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(opensemba_dgtd 
		maxwell-driver 
		mfem
		${GSL_LIBRARIES}
)

if (SEMBA_DGTD_ENABLE_OPENMP)
	message(STATUS "Linking opensemba_dgtd with OpenMP libraries")
	target_link_libraries(opensemba_dgtd OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()

message(STATUS "Creating build system for opensemba_rcs")

add_executable(opensemba_rcs "rcs_launcher.cpp")

target_link_libraries(opensemba_rcs
		maxwell-driver
		maxwell-components
		mfem
)

if (SEMBA_DGTD_ENABLE_OPENMP)
	message(STATUS "Linking opensemba_rcs with OpenMP libraries")
	target_link_libraries(opensemba_rcs OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()