cmake_minimum_required(VERSION 3.0)

message(STATUS "Building maxwell-evolution")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})
message(STATUS "MFEM_INCLUDE_DIRS=" ${MFEM_INCLUDE_DIRS})

add_library(maxwell-evolution STATIC 
	"Evolution3D.cpp"
	"EvolutionMethods.cpp"
)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(maxwell-evolution mfem mfemExtension Eigen3::Eigen)
