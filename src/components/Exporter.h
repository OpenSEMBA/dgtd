#pragma once

#include <mfem.hpp>

namespace maxwell {

using namespace mfem;

class NearToFarFieldExporter : public DataCollection 
{
public:
	NearToFarFieldExporter(std::string name);

};

}