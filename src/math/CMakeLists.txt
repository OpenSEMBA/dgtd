message(STATUS "Creating build system for maxwell-math")

add_library(maxwell-math STATIC 
	"Geometry.cpp"
)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(maxwell-math mfem)

if (SEMBA_DGTD_ENABLE_OPEN_MP)
	message(STATUS "Linking maxwell-math with OpenMP libraries")
	target_link_libraries(maxwell-math OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()
