message(STATUS "Creating build system for mfem_tests")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(FFTW3 CONFIG REQUIRED)

include_directories(
	./
)

add_executable(rcs_tests
	"RCSTest.cpp"
)

target_link_libraries(rcs_tests 
	mfem 
	maxwell-mfemExtension
	maxwell-components
	maxwell-driver
	Eigen3::Eigen
	GTest::gtest 
	GTest::gtest_main
	FFTW3::fftw3
)