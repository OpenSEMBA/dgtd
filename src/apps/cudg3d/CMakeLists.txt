cmake_minimum_required(VERSION 2.8.11)

project(cudg3d CXX)

add_subdirectory(../../core obj/core)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

add_executable(cudg3d cudg3d.cpp)

target_link_libraries(cudg3d cudg3d_core
							 pthread)