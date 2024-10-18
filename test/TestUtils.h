#pragma once

#include "math/PhysicalConstants.h"
#include <gtest/gtest.h>

static std::string testDataFolder()     { return "./testData/"; }
static std::string gmshMeshesFolder()   { return testDataFolder() + "gmshMeshes/"; }
static std::string mfemMeshesFolder()   { return testDataFolder() + "mfemMeshes/"; }
static std::string mfemMeshes1DFolder() { return testDataFolder() + "mfemMeshes/1D/"; }
static std::string mfemMeshes2DFolder() { return testDataFolder() + "mfemMeshes/2D/"; }
static std::string mfemMeshes3DFolder() { return testDataFolder() + "mfemMeshes/3D/"; }
static std::string smbInputsFolder()    { return testDataFolder() + "smbInputs/"; }
static std::string maxwellInputsFolder(){ return testDataFolder() + "maxwellInputs/"; }

static std::string smbCase(const std::string& caseName)
{
	return smbInputsFolder() + caseName + "/" + caseName + ".smb.json";
}

static std::string maxwellCase(const std::string& caseName)
{
	return maxwellInputsFolder() + caseName + "/" + caseName + ".json";
}

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

