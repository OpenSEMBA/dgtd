#pragma once

#include <mfem.hpp>

#include <components/Probes.h>
#include <solver/SolverOptions.h>

namespace maxwell {

using namespace mfem;

class RCSCollection : public DataCollection {
public:
	RCSCollection(const NearToFarFieldProbe&, const SolverOptions&);

	void Save();
private:


};

}