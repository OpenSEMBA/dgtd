message(STATUS "Creating build system for case_manager")

find_package(nlohmann_json CONFIG REQUIRED)
	
add_library(maxwell-driver "driver.cpp" "rcs_driver.cpp")

target_link_libraries(maxwell-driver
	maxwell-solver
	mfem
	nlohmann_json
	${GSL_LIBRARIES}
)

target_include_directories(maxwell-driver PRIVATE ${GSL_INCLUDE_DIRS})

if (SEMBA_DGTD_ENABLE_OPENMP)
	message(STATUS "Linking maxwell-driver with OpenMP libraries")
	target_link_libraries(maxwell-driver OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()