cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

add_subdirectory(pybind11)

pybind11_add_module(maxwell_solver Solver_binding.cpp)