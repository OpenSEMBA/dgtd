project(test)
cmake_minimum_required(VERSION 2.8)

find_package( PkgConfig REQUIRED)
pkg_check_modules( gtest REQUIRED gtest)

add_sources(. SRCS)
add_executable(test) 

target_link_libraries(test gtest)