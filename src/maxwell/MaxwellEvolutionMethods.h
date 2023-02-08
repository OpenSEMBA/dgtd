#include "Types.h"
#include "mfem.hpp"
#include "Model.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"
#include "mfemExtension/LinearIntegrators.h"
#include "mfemExtension/LinearForm_IBFI.hpp"

namespace maxwell {

using namespace mfem;
using FiniteElementIBFIOperator = std::unique_ptr<mfemExtension::BilinearFormIBFI>;
using FiniteElementOperator = std::unique_ptr<mfemExtension::BilinearForm>;
using FiniteElementVector = std::unique_ptr<mfemExtension::LinearFormIBFI>;

FiniteElementOperator buildByMult					(const BilinearForm& op1, const BilinearForm& op2                    , FiniteElementSpace&);
FiniteElementIBFIOperator buildIBFIByMult			(const BilinearForm& op1, const mfemExtension::BilinearFormIBFI& op2, FiniteElementSpace& fes);

FiniteElementOperator buildInverseMassMatrix		(const FieldType&, const Model&, FiniteElementSpace&);
FiniteElementOperator buildDerivativeOperator		(const Direction&, FiniteElementSpace&);
FiniteElementOperator buildFluxOperator				(const FieldType&, const std::vector<Direction>&, bool usePenaltyCoefficients, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
FiniteElementOperator buildFluxJumpOperator			(const FieldType&, const std::vector<Direction>&, bool usePenaltyCoefficients, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
FiniteElementOperator buildFluxOperator				(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
FiniteElementOperator buildPenaltyOperator			(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
FiniteElementIBFIOperator buildFunctionOperator		(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);

FiniteElementOperator buildFluxOperator1D				(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
FiniteElementOperator buildPenaltyOperator1D			(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
FiniteElementIBFIOperator buildFluxFunctionOperator1D	(Model&, FiniteElementSpace&);
FiniteElementIBFIOperator buildPenaltyFunctionOperator1D(Model&, FiniteElementSpace&);

FiniteElementVector   buildBoundaryFunctionVector1D		(Model&, FiniteElementSpace&);

TFSFOrientationCoefficient interiorBoundaryFaceCoefficient(const BdrCond&);

FluxCoefficient interiorFluxCoefficient();
FluxCoefficient interiorPenaltyFluxCoefficient		(const MaxwellEvolOptions&);
FluxCoefficient boundaryFluxCoefficient				(const FieldType&, const BdrCond&);
FluxCoefficient boundaryPenaltyFluxCoefficient		(const FieldType&, const BdrCond&, const MaxwellEvolOptions&);

FieldType altField(const FieldType& f);

}