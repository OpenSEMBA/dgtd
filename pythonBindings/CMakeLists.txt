cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	message(FATAL_ERROR "pybindings do not work in debug mode." )
endif()

add_subdirectory(pybind11)

pybind11_add_module(maxwell_solver SHARED maxwell_solver_binding.cpp)

target_link_libraries(maxwell_solver PRIVATE maxwell)