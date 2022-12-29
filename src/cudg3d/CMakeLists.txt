cmake_minimum_required(VERSION 3.0)

#add_subdirectory(communications)
add_subdirectory(dg)
add_subdirectory(integrator)

find_package(nlohmann_json CONFIG REQUIRED)
find_package(mfem CONFIG REQUIRED)

message(STATUS "Making build system for dgtd_solver.")
message(STATUS "mfem include directories: " ${MFEM_INCLUDE_DIRS})
include_directories(${MFEM_INCLUDE_DIRS})

add_library(dgtd_solver 
	"Cudg3d.cpp" 	
) 
target_link_libraries(dgtd_solver 
	dgtd_dg 
	dgtd_integrator 
	opensemba_core 
	nlohmann_json
)
