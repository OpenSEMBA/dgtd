message(STATUS "Creating build system for case_manager")

find_package(nlohmann_json CONFIG REQUIRED)

include_directories(./)

add_executable (cases_tests
    $<TARGET_OBJECTS:TestLauncher>
	"CasesTest.cpp"
	../FieldUtils.cpp)

if(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS)
	message(STATUS "Appending ExtensiveCasesTests to cases_tests")
	target_sources(cases_tests PRIVATE
		"ExtensiveCasesTest.cpp")
endif()

if(SEMBA_DGTD_ENABLE_L2_ANALYSIS_TESTS)
	message(STATUS "Appending AnalyticalStudyTests and L2ErrorAnalysis to cases_tests")
	target_sources(cases_tests PRIVATE
		"AnalyticalStudyTests.cpp"
		"L2ErrorAnalysis.cpp")
endif()

target_include_directories(cases_tests PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(cases_tests maxwell-driver GTest::gtest mfem ${GSL_LIBRARIES})