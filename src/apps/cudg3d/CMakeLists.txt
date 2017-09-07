CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)
PROJECT(cudg3d)

FIND_PACKAGE( PkgConfig REQUIRED)

# Moving the binary files to the directory bin/
SET (EXECUTABLE_OUTPUT_PATH ${Main_BINARY_DIR}/bin
     CACHE PATH "Directory for Executables")

OPTION(RUN_TESTS "Compile tests" ON)
if(RUN_TESTS)
	FIND_PACKAGE( PkgConfig REQUIRED)
	PKG_CHECK_MODULES( gtest REQUIRED gtest)
	add_subdirectory(src/apps/test)
endif()

OPTION(OPENSEMBA "Inlcude OpenSEMBA library" ON)
if(OPENSEMBA)
	#include_directories(${CMAKE_CURRENT_LIST_DIR}
        #                   ${CMAKE_CURRENT_LIST_DIR}/external/)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/OpenSEMBA/src ${CMAKE_CURRENT_BINARY_DIR}/opensemba)
endif()

ADD_EXECUTABLE(Main Main.cpp)
TARGET_LINK_LIBRARIES(Main ${gsl_LIBRARIES} pthread)
