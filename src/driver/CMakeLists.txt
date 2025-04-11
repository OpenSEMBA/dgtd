message(STATUS "Creating build system for case_manager")

find_package(nlohmann_json CONFIG REQUIRED)
	
add_library(maxwell-driver "driver.cpp")

target_link_libraries(maxwell-driver
	maxwell-solver
	mfem
	nlohmann_json
)

if (SEMBA_DGTD_ENABLE_OPEN_MP)
	message(STATUS "Linking maxwell-driver with OpenMP libraries")
	target_link_libraries(maxwell-driver OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()