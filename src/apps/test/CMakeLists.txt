cmake_minimum_required(VERSION 2.8)

add_sources(. SRCS)
add_executable(test) 

target_link_libraries(semba semba_semba semba_util)