#include "Probes.h"
#include <algorithm>

using namespace mfem;

namespace maxwell {

void checkPointsHaveSameSize(const Points& points)
{
	for (int i = 0; i < points.size() - 1; i++) {
		if (points.at(i).size() != points.at(i + 1).size()) {
			throw std::exception("Points do not have same dimension.");
		}
	}
}

void checkPointsAreNotEmpty(const std::vector<Point>& points) 
{
	for (int i = 0; i < points.size(); i++) {
		if (points.at(i).size() == 0) {
			throw std::exception("Empty subvector in points at position " + i);
		}
	}
}

PointsProbe::PointsProbe(const FieldType& ft, const Direction& d, const Points& points) :
	fieldToExtract_{ ft },
	directionToExtract_{ d },
	points_{points}
{
	if (points.size() == 0) {
		throw std::exception("Empty points vector.");
	}

	checkPointsAreNotEmpty(points);
	checkPointsHaveSameSize(points);
}

}