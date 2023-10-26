#pragma once

#include "TestUtils.h"
#include "components/probes.h"

namespace maxwell {

static Probes buildProbesWithAnExportProbe(int visSteps = 1) 
{
#ifdef MAXWELL_TEST_ALLOW_PARAVIEW_EXPORT
	Probes r{ {}, { ExporterProbe{ getTestCaseName()} } };
	r.exporterProbes[0].visSteps = visSteps;
	return r;
#else
	return {{},{}};
#endif
}

}