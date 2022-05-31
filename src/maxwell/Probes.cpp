#include "Probes.h"

using namespace mfem;

namespace maxwell {
Probe::Probe(const FieldType& ft, const Direction& d, std::vector<std::vector<double>>& points) :
	fieldToExtract_(ft),
	directionToExtract_(d)
{

	integPointMat_.SetSize(points.at(0).size(), points.size());
	for (int i = 0; i < points.at(i).size(); i++) {
		for (int j = 0; j < points.size(); j++) {
			integPointMat_.Elem(i, j) = points.at(j).at(i);
		}
	}
}

}