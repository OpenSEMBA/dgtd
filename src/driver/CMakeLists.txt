cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for case_manager")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(nlohmann_json CONFIG REQUIRED)

add_library(maxwell-driver "driver.cpp"
)

target_link_libraries(maxwell-driver
	mfem
	maxwell-solver
	nlohmann_json)