#pragma once

#include <mfem.hpp>

#include "components/Model.h"
#include "EvolutionOptions.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"
#include "mfemExtension/LinearIntegrators.h"
#include "mfemExtension/LinearForm_IBFI.hpp"

#include "components/DGOperatorFactory.h"

#include "math/EigenMfemTools.h"
#include "solver/SourcesManager.h"

namespace maxwell {

	using namespace mfem;

	const FieldGridFuncs evalTimeVarFunction(const Time time, SourcesManager& sm);

	std::vector<int> calcOffsetCoeff(const std::vector<FieldType>&, const std::vector<Direction>&);

	void allocateDenseInEigen(DenseMatrix* bilForm, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const std::vector<Direction> d, const double sign = 1.0);

}
