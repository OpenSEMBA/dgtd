message(STATUS "Creating build system for maxwell-components")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(nlohmann_json CONFIG REQUIRED)

add_library(maxwell-components STATIC 
	"Material.cpp"
	"Model.cpp" 
	"Sources.cpp"
	"EigenvalueEstimator.cpp"
	"OptimizationManager.cpp"
	"SubMesher.cpp"
	"RCSManager.cpp" 
	"RCSExporter.cpp"
)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(maxwell-components mfem Eigen3::Eigen nlohmann_json)
