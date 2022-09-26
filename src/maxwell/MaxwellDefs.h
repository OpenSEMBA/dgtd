#include "Types.h"
#include "SolverOptions.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;
using FiniteElementOperator = std::unique_ptr<BilinearForm>;

Vector buildNVector(const Direction& d, const FiniteElementSpace& fes);

FiniteElementOperator buildByMult(const BilinearForm& op1,const BilinearForm& op2, FiniteElementSpace& fes);
FiniteElementOperator buildInverseMassMatrix(const FieldType& f, const Model& model, FiniteElementSpace& fes);
FiniteElementOperator buildDerivativeOperator(const Direction& d, FiniteElementSpace& fes);
FiniteElementOperator buildFluxOperator(const FieldType& f, const Direction& d, bool usePenaltyCoefficients, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);

FluxCoefficient interiorFluxCoefficient();
FluxCoefficient interiorPenaltyFluxCoefficient(const MaxwellEvolOptions& opts);
FluxCoefficient boundaryFluxCoefficient(const FieldType& f, const BdrCond& bdrC);
FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType& f, const BdrCond& bdrC, const MaxwellEvolOptions& opts);

FieldType altField(const FieldType& f);
}