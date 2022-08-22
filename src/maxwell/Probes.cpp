#include "Probes.h"
#include <algorithm>

using namespace mfem;

namespace maxwell {

ExporterProbe::ExporterProbe(int order) :
	order_{ order }
{}



const bool PointsProbe::verifyEntryVectorsSameSize(std::vector<std::vector<double>>& points) const
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

const void PointsProbe::verifyEntrySubvectorsNotEmpty(std::vector<std::vector<double>>& points) const
{
	for (int i = 0; i < points.size(); i++) {
		if (points.at(i).size() == 0) {
			throw std::exception("Empty subvector in points at position " + i);
		}
	}
}

const void PointsProbe::buildIntegPointMat(std::vector<std::vector<double>>& points)
{
	integPointMat_.SetSize((int) points.at(0).size(),(int) points.size());
	if ((int) points.at(0).size() == 1 && (int) points.size() == 1) {
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

PointsProbe::PointsProbe(const FieldType& ft, const Direction& d, std::vector<std::vector<double>>& points) :
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