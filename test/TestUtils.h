#pragma once

#include <gtest/gtest.h>

static std::string testDataFolder()   { return "./testData/"; }
static std::string gmshMeshesFolder() { return testDataFolder() + "gmshMeshes/"; }
static std::string mfemMeshesFolder() { return testDataFolder() + "mfemMeshes/"; }
static std::string smbInputsFolder()  { return testDataFolder() + "smbInputs/"; }

static std::string getCaseName()
{
	return ::testing::UnitTest::GetInstance()->current_test_info()->name();
}

static std::string getTestCaseName() 
{
	std::string suiteName{
		::testing::UnitTest::GetInstance()->current_test_info()->test_suite_name()
	};	
	return suiteName + "." + getCaseName();
}

