message(STATUS "Creating build system for case_manager")

find_package(nlohmann_json CONFIG REQUIRED)

add_library(maxwell-driver "driver.cpp"
)

target_link_libraries(maxwell-driver
	maxwell-solver
	mfem
	nlohmann_json
)