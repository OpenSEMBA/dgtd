cmake_minimum_required(VERSION 3.0)

add_library(dgtd_dg 
	"Evolution.cpp"
	"VolumeModel.cpp"
 "CellTet.h" "Field.h")

target_link_libraries(dgtd_dg opensemba_core)
