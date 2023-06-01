cmake_minimum_required(VERSION 3.0)

set(LIBRARY_NAME "maxwell-adapter")

message(STATUS "Creating build system for " ${LIBRARY_NAME})

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

set(OPENSEMBA_ROOT_DIR ${PROJECT_SOURCE_DIR}/external/opensemba-core/)
add_subdirectory(${OPENSEMBA_ROOT_DIR} opensemba-core)

add_library(${LIBRARY_NAME} STATIC 
	"OpensembaAdapter.cpp"
)
include_directories(${OPENSEMBA_ROOT_DIR}/src/)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(${LIBRARY_NAME} mfem opensemba) 
