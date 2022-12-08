cmake_minimum_required(VERSION 3.0)

#add_subdirectory(jacobi)
#add_subdirectory(boundaryCondition)


find_package(nlohmann_json CONFIG REQUIRED)

add_library(dgtd_solver 
	"Cudg3d.cpp" 	
)

target_link_libraries(dgtd_solver nlohmann_json)
