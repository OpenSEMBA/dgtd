#include "Probes.h"

using namespace mfem;

namespace maxwell {
Probe::Probe(const FieldType& ft, const Direction& d, const DenseMatrix& mat) :
	fieldToExtract_(ft),
	directionToExtract_(d),
	integPointMat_(mat)
{}

}