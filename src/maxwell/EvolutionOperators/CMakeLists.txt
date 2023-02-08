cmake_minimum_required(VERSION 3.0)

project(maxwell)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(EvolutionOP STATIC 
"MaxwellEvolution1D_FA.cpp"
)

target_link_libraries(EvolutionOP mfem Eigen3::Eigen)
