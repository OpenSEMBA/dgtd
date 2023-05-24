#pragma once

#include "Utils.h"
#include "components/probes.h"

namespace maxwell {

static Probes buildProbesWithAnExportProbe() {
		return { {}, { ExporterProbe{ getTestCaseName()} } };
}

}