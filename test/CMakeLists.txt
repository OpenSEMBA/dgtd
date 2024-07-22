cmake_minimum_required(VERSION 3.5)

include_directories(
	./
	${maxwell_INCLUDE_DIRS}
)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_subdirectory(mfem)
add_subdirectory(hesthavenComparison)

add_subdirectory(maxwell)

add_subdirectory(cases)
add_subdirectory(rcs)

add_test(mfem ${CMAKE_BINARY_DIR}/bin/mfem_tests)
add_test(hesthavenComparison ${CMAKE_BINARY_DIR}/bin/hesthavenComparison_tests)
add_test(maxwell ${CMAKE_BINARY_DIR}/bin/maxwell_tests)

add_test(rcs ${CMAKE_BINARY_DIR}/bin/rcs_tests)
add_test(maxwell_solver ${CMAKE_BINARY_DIR}/bin/maxwell_solver_tests)
add_test(cases ${CMAKE_BINARY_DIR}/bin/cases_tests)
