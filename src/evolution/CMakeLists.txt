message(STATUS "Creating build system for maxwell-evolution")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(maxwell-evolution STATIC 
	"MaxwellEvolution.cpp"
	"MaxwellEvolutionMethods.cpp"
	"Fields.cpp" 
	"HesthavenEvolutionMethods.cpp"
	"HesthavenEvolution.cpp" 
	"GlobalEvolution.cpp")

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(maxwell-evolution maxwell-components mfem maxwell-mfemExtension Eigen3::Eigen)

if (SEMBA_DGTD_ENABLE_OPEN_MP)
	message(STATUS "Linking maxwell-evolution with OpenMP libraries")
	target_link_libraries(maxwell-evolution OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()
