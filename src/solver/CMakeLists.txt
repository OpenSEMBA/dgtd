message(STATUS "Creating build system for maxwell-solver")

if(SEMBA_DGTD_ENABLE_TIMER_INFORMATION)
	add_definitions(-DSHOW_TIMER_INFORMATION)
endif()

add_library(maxwell-solver STATIC 
	"Solver.cpp" 
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
	target_link_libraries(maxwell-solver OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()

if(SEMBA_DGTD_ENABLE_CUDA)
	add_definitions(-DSEMBA_DGTD_ENABLE_CUDA)
endif()
