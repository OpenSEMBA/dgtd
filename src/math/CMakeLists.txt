cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for maxwell-math")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell-math STATIC 
	"Geometry.cpp"
)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(maxwell-math mfem)
