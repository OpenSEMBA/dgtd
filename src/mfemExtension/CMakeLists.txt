cmake_minimum_required(VERSION 3.0)

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(mfemExtension STATIC 
	"BilinearForm_IBFI.cpp"
	"LinearForm_IBFI.cpp"
	"BilinearIntegrators.cpp"
	"LinearIntegrators.cpp" 
	"IntegratorFunctions.cpp" 
)

target_link_libraries(mfemExtension mfem)
