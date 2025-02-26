#include "MaxwellEvolutionMethods.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

const FieldGridFuncs evalTimeVarFunction(const Time time, SourcesManager& sm)
{
	auto res{ sm.evalTimeVarField(time, sm.getGlobalTFSFSpace()) };
	auto func_g_sf = res;
	sm.markDoFSforTForSF(res, true);
	{
		if (sm.getTFSFSubMesher().getSFSubMesh() != NULL) {
			sm.markDoFSforTForSF(func_g_sf, false);
			for (int f : {E, H}) {
				for (int x{ 0 }; x <= Z; x++) {
					res[f][x] -= func_g_sf[f][x];
					res[f][x] *= 0.5;
				}
			}
		}
	}
	return res;
}

std::vector<int> calcOffsetCoeff(const std::vector<FieldType>& f, const std::vector<Direction>& d)
{
	std::vector<int> res(2);
	if (d.size() == 1) {
		if (f[0] == E) {
			res[0] = d[0];
			res[1] = d[0];
		}
		else {
			res[0] = 3 + d[0];
			res[1] = 3 + d[0];
		}
	}
	else if (f[0] == f[1]) {
		if (f[0] == E) {
			res[0] = d[0];
			res[1] = d[1];
		}
		else {
			res[0] = 3 + d[0];
			res[1] = 3 + d[1];
		}
	}
	else if (f[0] != f[1]) {
		if (f[0] == E) {
			res[0] = d[0];
			res[1] = 3 + d[1];
		}
		else {
			res[0] = 3 + d[0];
			res[1] = d[1];
		}
	}
	else {
		throw std::runtime_error("Wrong input in method, check direction or field type vectors.");
	}
	return res;
}

void allocateDenseInEigen(DenseMatrix* bilMat, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const std::vector<Direction> d, const double sign)
{
	auto offset = bilMat->Height();
	auto offsetCoeff{ calcOffsetCoeff(f,d) };

	for (int i = 0; i < bilMat->Height(); ++i) {
		for (int j = 0; j < bilMat->Width(); ++j) {
			if (bilMat->Elem(i, j) != 0.0) {
				res.coeffRef(i + offset * offsetCoeff[0], j + offset * offsetCoeff[1]) += sign * bilMat->Elem(i, j);
			}
		}
	}
}

}