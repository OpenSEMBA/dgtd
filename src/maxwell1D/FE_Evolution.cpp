#include "FE_Evolution.h"

namespace Maxwell1D {

FE_Evolution::FE_Evolution(FiniteElementSpace* fes) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	fes_(fes),
	MInv_(buildInverseMassMatrix()),
	KxE_(buildDerivativeAndFluxOperator(X, Electric)),
	KxH_(buildDerivativeAndFluxOperator(X, Magnetic))
{}

std::unique_ptr<BilinearForm> FE_Evolution::buildInverseMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}

std::unique_ptr<BilinearForm> FE_Evolution::buildDerivativeAndFluxOperator(
	const Direction& d, const FieldType& ft) const
{
	assert(d == X, "Incorrect argument for direction.");

	auto K = std::make_unique<BilinearForm>(fes_);

	ConstantCoefficient one(1.0);

	K->AddDomainIntegrator(
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
		K->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		K->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}
	else
	{
		alpha = -1.0;
		beta = 0.0;
		K->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		K->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}

	K->Assemble();
	K->Finalize();

	return K;
}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	Vector eOld(x.GetData(),                    fes_->GetNDofs());
	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());

	GridFunction eNew(fes_, &y[0]);
	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);

	Vector aux(MInv_->Height());

	// Update E.
	KxH_->Mult(hOld, aux);
	MInv_->Mult(aux, eNew);
	 
	// Update H.
	KxE_->Mult(eOld, aux);
	MInv_->Mult(aux, hNew);
}

}