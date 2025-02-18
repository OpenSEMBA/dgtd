message(STATUS "Creating build system for maxwell-mfemExtension")

add_library(maxwell-mfemExtension STATIC 
	#"BilinearForm_IBFI.cpp"
	#"LinearForm_IBFI.cpp"
	"BilinearIntegrators.cpp"
	"LinearIntegrators.cpp" 
	"IntegratorFunctions.cpp" 
)

target_link_libraries(maxwell-mfemExtension mfem)
