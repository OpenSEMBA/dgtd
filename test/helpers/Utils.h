#pragma once

#include <gtest/gtest.h>
#include "components/Probes.h"

namespace maxwell {

std::string getTestCaseName();

Probes buildProbesWithAnExportProbe();


}