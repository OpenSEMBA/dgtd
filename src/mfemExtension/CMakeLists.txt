cmake_minimum_required(VERSION 3.5)

message(STATUS "Creating build system for maxwell-mfemExtension")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(maxwell-mfemExtension STATIC 
	#"BilinearForm_IBFI.cpp"
	#"LinearForm_IBFI.cpp"
	"BilinearIntegrators.cpp"
	"LinearIntegrators.cpp" 
	"IntegratorFunctions.cpp" 
)

target_link_libraries(maxwell-mfemExtension mfem)
