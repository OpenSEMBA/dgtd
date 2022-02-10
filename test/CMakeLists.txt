cmake_minimum_required(VERSION 3.0)

add_sources(. SRCS)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(mfem CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_executable(maxwell_tests ${SRCS})
target_link_libraries(maxwell_tests maxwell mfem GTest::gtest GTest::gtest_main) 


