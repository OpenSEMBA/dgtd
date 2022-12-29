cmake_minimum_required(VERSION 3.0)

add_library(dgtd_integrator 
	"TimeIntegrator.cpp"
	"LSERK4.cpp" 	
	#"LF2.cpp"
)

target_link_libraries(dgtd_integrator opensemba_core)
