#include "Types.h"
#include "mfem.hpp"
#include "Model.h"
#include "mfemExtension/BilinearIntegrators.h"


namespace maxwell {

using namespace mfem;
using FiniteElementOperator = std::unique_ptr<BilinearForm>;

Vector buildNVector(const Direction& d, const FiniteElementSpace& fes);

FiniteElementOperator buildFluxOperator1D(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);
FiniteElementOperator buildPenaltyOperator1D(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);

FluxCoefficient interiorCenteredFluxCoefficient1D();
FluxCoefficient interiorPenaltyFluxCoefficient1D(const MaxwellEvolOptions& opts);
FluxCoefficient boundaryCenteredFluxCoefficient1D(const FieldType& f, const BdrCond& bdrC);
FluxCoefficient boundaryPenaltyFluxCoefficient1D(const FieldType& f, const BdrCond& bdrC, const MaxwellEvolOptions& opts);

}