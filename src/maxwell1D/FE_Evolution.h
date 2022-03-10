#pragma once

#include "mfem.hpp"

namespace Maxwell1D {

using namespace mfem;

class FE_Evolution : public TimeDependentOperator
{
public:

	FE_Evolution(FiniteElementSpace* fes);
	virtual void Mult(Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:

	typedef std::size_t Direction;
	typedef std::size_t FieldType;

	const Direction X = 0;

	const FieldType Electric = 0;
	const FieldType Magnetic = 1;

	FiniteElementSpace* fes_;
	Vector z_;

	BilinearForm& MInv_;
	BilinearForm& KxE_;
	BilinearForm& KxH_;

	BilinearForm& buildInverseMassMatrix() const;
	BilinearForm& buildDerivativeAndFluxOperator(
		const Direction& d, const FieldType& ft) const;

};

FE_Evolution::FE_Evolution(FiniteElementSpace* fes) :
	fes_(fes),
	MInv_(buildInverseMassMatrix()),
	KxE_(buildDerivativeAndFluxOperator(X, Electric)),
	KxH_(buildDerivativeAndFluxOperator(X, Magnetic)),
	TimeDependentOperator(MInv_.Height()),
	z_(MInv_.Height())
{}

BilinearForm& FE_Evolution::buildInverseMassMatrix() const
{
	BilinearForm MInv = BilinearForm(fes_);
	MInv.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	MInv.Assemble();
	MInv.Finalize();
	return MInv;
}

BilinearForm& FE_Evolution::buildDerivativeAndFluxOperator(
	const Direction& d, const FieldType& ft) const
{
	assert(d == X, "Incorrect argument for direction.");

	auto K = BilinearForm(fes_);

	ConstantCoefficient one(1.0);

	K.AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, d)));

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	double alpha;
	double beta;

	if (ft == Electric)
	{
		alpha = -1.0;
		beta = 0.0;
		K.AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		K.AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}
	else
	{
		alpha = -1.0;
		beta = 0.0;
		K.AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		K.AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}

	K.Assemble();
	K.Finalize();

	return K;
}

void FE_Evolution::Mult(Vector& x, Vector& y) const
{
	GridFunction eOld(fes_,&x[0]);
	GridFunction hOld(fes_,&x[fes_->GetNDofs()]);

	GridFunction eNew(fes_,&y[0]);
	GridFunction hNew(fes_,&y[fes_->GetNDofs()]);

	Vector aux(MInv_.Height());

	// Update E.
	KxH_.Mult(hOld, aux);
	MInv_.Mult(aux, eNew);
	 
	// Update H.
	KxE_.Mult(eOld, aux);
	MInv_.Mult(aux, hNew);
}

}