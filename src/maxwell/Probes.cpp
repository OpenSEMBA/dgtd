#include "Probes.h"

using namespace mfem;

namespace maxwell {
Probe::Probe(const FieldType& ft, const Direction& d, std::vector<std::vector<double>>& points) :
	fieldToExtract_(ft),
	directionToExtract_(d)
{
	integPointMat_.SetSize(points.size());
	for (int i = 0; i < points.size(); i++) {
		for (int j = 0; j < points.at(i).size(); j++) {
			integPointMat_.Elem(i, j) = points.at(i).at(j);
		}
	}
}

}