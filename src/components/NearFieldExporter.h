#pragma once

#include <mfem.hpp>

#include <components/Probes.h>
#include <solver/SolverOptions.h>

namespace maxwell {

using namespace mfem;

class NearFieldCollection : public DataCollection {
public:
	NearFieldCollection(const SurfaceProbe&, const SolverOptions&);

	void Save();
private:


};

}