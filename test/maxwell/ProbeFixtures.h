#pragma once

#include "TestUtils.h"
#include "components/Probes.h"

namespace maxwell {

static Probes buildProbesWithAnExportProbe(int visSteps = 1) 
{
	Probes r{ {}, { ExporterProbe{ getTestCaseName()} } };
	r.exporterProbes[0].visSteps = visSteps;
	return r;
}

static Probes buildProbesEmpty()
{
	return { {}, {} };
}

}
