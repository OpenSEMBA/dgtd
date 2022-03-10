#pragma once

#include "mfem.hpp"

namespace Maxwell1D {

using namespace mfem;

class FE_Evolution : public TimeDependentOperator
{
public:
	FE_Evolution(
		BilinearForm& invM, 
		BilinearForm& Ke, 
		BilinearForm& Kh);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:
	BilinearForm &invM_, &Ke_, &Kh_;
	Vector z_;
};

FE_Evolution::FE_Evolution(BilinearForm& invM, BilinearForm& Ke, BilinearForm& Kh): 
	TimeDependentOperator(2*invM.Height()), 
	invM_(invM), 
	Ke_(Ke), 
	Kh_(Kh),
	z_(2*invM.Height())
{}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	// Update E.
	KxH_->Mult(Hy_, aux);
	MInv_->Mult(aux, ezNew);
	
	// Update H.
	KxE_->Mult(Ez_, aux);
	MInv_->Mult(aux, hyNew);

	Ez_ = ezNew;
	Hy_ = hyNew;
}

}