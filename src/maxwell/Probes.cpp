#include "Probes.h"
#include <algorithm>

using namespace mfem;

namespace maxwell {

const bool Probe::verifyEntryVectorsSameSize(std::vector<std::vector<double>>& points) const
{
	bool res{};

	if (points.size() == 1) {
		return true;
	}

	for (int i = 0; i < points.size() - 1; i++) {
		res = points.at(i).size() == points.at(i + 1).size();
		if (res == false)
			break;
	}
	return res;
}

const void Probe::verifyEntrySubvectorsNotEmpty(std::vector<std::vector<double>>& points) const
{
	for (int i = 0; i < points.size(); i++) {
		if (points.at(i).size() == 0) {
			throw std::exception("Empty subvector in points at position " + i);
		}
	}
}

const void Probe::buildIntegPointMat(std::vector<std::vector<double>>& points)
{
	integPointMat_.SetSize(points.at(0).size(), points.size());
	if (points.at(0).size() == 1 && points.size() == 1) {
		integPointMat_.Elem(0, 0) = points.at(0).at(0);
	}
	else {
		for (int i = 0; i < points.at(i).size(); i++) {
			for (int j = 0; j < points.size(); j++) {
				integPointMat_.Elem(i, j) = points.at(j).at(i);
			}
		}
	}
}

Probe::Probe(const FieldType& ft, const Direction& d, std::vector<std::vector<double>>& points) :
	fieldToExtract_(ft),
	directionToExtract_(d)
{
	if (points.size() == 0) {
		throw std::exception("Empty points vector.");
	}

	verifyEntrySubvectorsNotEmpty(points);

	bool ssize = verifyEntryVectorsSameSize(points);

	if (ssize) {
		buildIntegPointMat(points);
	}
	else {
		throw std::exception("Vectors in points do not have the same dimensions.");
	}
}
}
