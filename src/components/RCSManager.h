#pragma once

#include <mfem.hpp>
#include <components/Types.h>
#include <iostream>
#include <components/Probes.h>


namespace maxwell {

using namespace mfem;

using FieldToGridFunction = std::map<std::string, GridFunction>;
using FieldsToTime = std::pair<FieldToGridFunction, Time>;

class RCSManager {
public:


	RCSManager(const std::string& path, const NearToFarFieldProbe&);

	void update(FieldsToTime&);

private:

	void initFieldsRCS(const std::string&);
	void calculateRCS(FieldsToTime&);
	GridFunction getGridFunction(const std::string& path, const FieldType&, const Direction&);

	std::unique_ptr<Mesh> m_;
	FieldToGridFunction fieldsRCS_;
	std::string basePath_;

};

}
