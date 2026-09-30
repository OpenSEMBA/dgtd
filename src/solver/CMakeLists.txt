message(STATUS "Creating build system for maxwell-solver")

find_package(Eigen3 3.3 REQUIRED CONFIG)

add_library(maxwell-solver STATIC 
	"Solver.cpp"
	"SolverExtension.cpp"
	"ProbesManager.cpp" 
	"SourcesManager.cpp")

target_link_libraries(maxwell-solver 
	maxwell-evolution 
	maxwell-components 
	maxwell-math
	mfem
)

if (SEMBA_DGTD_ENABLE_OPENMP)
	message(STATUS "Linking maxwell-solver with OpenMP libraries")
	target_link_libraries(maxwell-solver Eigen3::Eigen OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()

if(SEMBA_DGTD_ENABLE_CUDA)
	add_definitions(-DSEMBA_DGTD_ENABLE_CUDA)
	target_sources(maxwell-solver PRIVATE "SourcesManager.cu")
	set_source_files_properties(SourcesManager.cu PROPERTIES LANGUAGE CUDA)
	set_target_properties(maxwell-solver PROPERTIES
		CUDA_STANDARD 17
		CUDA_STANDARD_REQUIRED ON
	)
endif()
