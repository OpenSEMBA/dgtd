#include <mfem.hpp>
#include <filesystem>
#include <fstream>

#include "RCSManager.h"
#include "SubMesher.h"
#include "../test/TestUtils.h"
#include "math/PhysicalConstants.h"

namespace maxwell {

using namespace mfem;

class RCSDataExtractor {
public:
	RCSDataExtractor(const std::string data_folder, const std::string case_name);
};

}