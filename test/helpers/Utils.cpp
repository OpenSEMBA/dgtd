#include "Utils.h"

namespace maxwell {

Probes buildProbesWithAnExportProbe()
{
	return { {}, { ExporterProbe{getTestCaseName()} } };
}

std::string getTestCaseName()
{
	std::string caseName{
		::testing::UnitTest::GetInstance()->current_test_info()->test_suite_name()
	};
	std::string name{
		::testing::UnitTest::GetInstance()->current_test_info()->name()
	};
	return caseName + "." + name;
}

}