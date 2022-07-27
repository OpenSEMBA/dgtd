#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"
#include "Types.h"
#include "Model.h"
#include "Sources.h"

#include <array>

namespace maxwell {

class FiniteElementEvolution : public TimeDependentOperator {
public:

	typedef std::unique_ptr<BilinearForm> FiniteElementOperator;

	struct Options {
		FluxType fluxType = FluxType::Upwind;
		DisForm disForm = DisForm::Weak;
	};

	static const std::size_t numberOfFieldComponents = 2;
	static const std::size_t numberOfMaxDimensions = 3;

	FiniteElementEvolution(FiniteElementSpace* fes, Options options, Model& model, Sources& sources);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FiniteElementEvolution() = default;

	const FiniteElementOperator& getInvMassStiffness(
		const FieldType& f, const Direction& dir) const 
	{ return MS_[f][dir]; }
	const FiniteElementOperator& getInvMassFlux(
		const FieldType& f1, const FieldType& f2, const Direction& dir) const 
	{ return MF_[f1][f2][dir]; }
	const FiniteElementOperator& getInvMassPenalty(
		const FieldType& f1, const FieldType& f2, const Direction& dir) const 
	{ return MP_[f1][f2][dir]; }
	const FiniteElementOperator& getInvMassNoDirFlux(
		const FieldType& f1, const FieldType& f2) const
	{ return MNND_[f1][f2]; }
	const FiniteElementOperator& getInvMassOneDirFlux(
		const FieldType& f1, const FieldType& f2, const Direction& dir) const 
	{ return MNOD_[f1][f2][dir]; }
	const FiniteElementOperator& getInvMassTwoDirFlux(
		const FieldType& f1, const FieldType& f2, const Direction& dir1, const Direction& dir2) const 
	{ return MNTD_[f1][f2][dir1][dir2]; }

private:
	struct FluxCoefficient {
		double alpha;
		double beta;
	};

	FiniteElementSpace* fes_;
	Options opts_;
	Model model_;
	Sources sources_;

	std::array<std::array<FiniteElementOperator, 3>, 2> MS_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MF_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MP_;
	std::array<std::array<FiniteElementOperator, 2>, 2> MNND_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MNOD_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MNTD_;

	Vector buildNVector(const Direction& d) const;
	Vector buildPieceWiseArgVector(const FieldType& f) const;
	
	FiniteElementOperator buildDerivativeOperator(const Direction&) const;
	FiniteElementOperator buildInverseMassMatrix(const FieldType&) const;
	FiniteElementOperator buildFluxOperator(const DisForm&, const FieldType&, const Direction&) const;
	FiniteElementOperator buildPenaltyOperator(const DisForm&, const FieldType&, const Direction&) const;
	FiniteElementOperator buildNormalFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms) const;

	FiniteElementOperator buildByMult(const BilinearForm*, const BilinearForm*) const;

	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorPenaltyFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const FieldType&, const BdrCond& bdrC) const;
	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&, const BdrCond& bdrC) const;
};

}