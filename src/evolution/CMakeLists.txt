cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for maxwell-evolution")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(maxwell-evolution STATIC 
	"Evolution.cpp"
	"EvolutionMethods.cpp"
	"Fields.cpp")

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(maxwell-evolution mfem maxwell-mfemExtension Eigen3::Eigen)
