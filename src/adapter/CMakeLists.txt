cmake_minimum_required(VERSION 3.0)

set(LIBRARY_NAME "maxwell-components")

message(STATUS "Creating build system for " ${LIBRARY_NAME})

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(${LIBRARY_NAME} STATIC 
	"OpensembaAdapter.cpp"
)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(${LIBRARY_NAME} mfem)
