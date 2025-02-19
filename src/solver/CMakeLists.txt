message(STATUS "Creating build system for maxwell-solver")

if(SEMBA_DGTD_ENABLE_TIMER_INFORMATION)
	add_definitions(-DSHOW_TIMER_INFORMATION)
endif()

add_library(maxwell-solver STATIC 
	"Solver.cpp" 
    "ProbesManager.cpp" 
	"SourcesManager.cpp" 
)

target_link_libraries(maxwell-solver 
	maxwell-evolution 
	maxwell-components 
	maxwell-math
	mfem
)
