message(STATUS "Creating build system for maxwell-mfemExtension")

add_library(maxwell-mfemExtension STATIC 
	#"BilinearForm_IBFI.cpp"
	#"LinearForm_IBFI.cpp"
	"BilinearIntegrators.cpp"
	"LinearIntegrators.cpp" 
	"IntegratorFunctions.cpp" 
)

target_link_libraries(maxwell-mfemExtension mfem)

if (SEMBA_DGTD_ENABLE_OPEN_MP)
	message(STATUS "Linking maxwell-mfemExtension with OpenMP libraries")
	target_link_libraries(maxwell-mfemExtension OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()
