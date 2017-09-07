cmake_minimum_required(VERSION 2.8.11)

project(cudg3d_core CXX)

include_directories(${CMAKE_CURRENT_LIST_DIR}
                    ${CMAKE_CURRENT_LIST_DIR}/../external/
                    ${OPENSEMBA_PATH}/src
                    ${OPENSEMBA_PATH}/src/core
                    ${OPENSEMBA_PATH}/external)

add_sources(. SRCS)

add_library(cudg3d_core STATIC ${SRCS})

