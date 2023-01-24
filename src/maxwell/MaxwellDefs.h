#include "Types.h"
#include "mfem.hpp"
#include "Model.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"

namespace maxwell {

using namespace mfem;
using FiniteElementOperator = std::unique_ptr<BilinearForm>;

FiniteElementOperator buildByMult(const BilinearForm& op1,const BilinearForm& op2, FiniteElementSpace& fes);
FiniteElementOperator buildInverseMassMatrix(const FieldType& f, const Model& model, FiniteElementSpace& fes);
FiniteElementOperator buildDerivativeOperator(const Direction& d, FiniteElementSpace& fes);
FiniteElementOperator buildFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms, bool usePenaltyCoefficients, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);
FiniteElementOperator buildFluxJumpOperator(const FieldType& f, const std::vector<Direction>& dirTerms, bool usePenaltyCoefficients, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);
FiniteElementOperator buildFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes);
FiniteElementOperator buildPenaltyOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);
FiniteElementOperator buildFunctionOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);


FluxCoefficient interiorFluxCoefficient();
FluxCoefficient interiorPenaltyFluxCoefficient(const MaxwellEvolOptions& opts);
FluxCoefficient boundaryFluxCoefficient(const FieldType& f, const BdrCond& bdrC);
FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType& f, const BdrCond& bdrC, const MaxwellEvolOptions& opts);

FieldType altField(const FieldType& f);

}