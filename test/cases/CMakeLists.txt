message(STATUS "Creating build system for case_manager")

find_package(nlohmann_json CONFIG REQUIRED)

include_directories(./)

add_executable (cases_tests
    $<TARGET_OBJECTS:TestLauncher>
	"CasesTest.cpp")

if(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS)
	message(STATUS "Appending ExtensiveCasesTests to cases_tests")
	target_sources(cases_tests PRIVATE
		"ExtensiveCasesTest.cpp")
endif()

target_link_libraries(cases_tests maxwell-driver GTest::gtest)