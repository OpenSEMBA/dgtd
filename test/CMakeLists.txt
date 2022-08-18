cmake_minimum_required(VERSION 3.0)

include_directories(
	./helpers
	${maxwell_solvers_INCLUDE_DIRS}
)

add_subdirectory(helpers)
add_subdirectory(maxwell)
add_subdirectory(mfem)
add_subdirectory(mfemHesthavenComparison)

