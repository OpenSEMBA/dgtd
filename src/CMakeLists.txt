cmake_minimum_required(VERSION 3.0)

add_sources(. SRCS)

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell STATIC ${SRCS})

