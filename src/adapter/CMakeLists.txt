cmake_minimum_required(VERSION 3.0)

set(LIBRARY_NAME "maxwell-adapter")

message(STATUS "Creating build system for " ${LIBRARY_NAME})

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(opensemba CONFIG REQUIRED)
message(STATUS "OPENSEMBA_FOUND=" ${OPENSEMBA_FOUND})
message(STATUS "OPENSEMBA_INCLUDE_DIR=" ${OPENSEMBA_INCLUDE_DIR})

add_library(${LIBRARY_NAME} STATIC 
	"OpensembaAdapter.cpp"
)
include_directories(${OPENSEMBA_INCLUDE_DIRS})

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(${LIBRARY_NAME} mfem)
