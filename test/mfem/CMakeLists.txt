message(STATUS "Creating build system for mfem_tests")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(FFTW3 CONFIG REQUIRED)

find_package(FFTW3f CONFIG REQUIRED)

find_package(FFTW3l CONFIG REQUIRED)

include_directories(
	./
)

add_executable(mfem_tests 
	"FiniteElementSpaceTest.cpp"
	"MeshTest.cpp" 
	"GridFunctionTest.cpp"
	"AlgebraTest.cpp"
	"Forms.cpp")

target_link_libraries(mfem_tests PRIVATE
	mfem 
	maxwell-mfemExtension
	maxwell-components 
	maxwell-driver
	Eigen3::Eigen
	FFTW3::fftw3
	FFTW3::fftw3f
	FFTW3::fftw3l
	GTest::gtest 
	GTest::gtest_main
)


