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
	mutable Vector e_, h_;

};

FE_Evolution::FE_Evolution(BilinearForm& invM, BilinearForm& Ke, BilinearForm& Kh): 
	TimeDependentOperator(invM.Height()), 
	invM_(invM), 
	Ke_(Ke), 
	Kh_(Kh),
	e_(invM.Height()), 
	h_(invM.Height())
{}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	// Update E.
	KxH_->Mult(Hy_, aux);
	MInv_->Mult(aux, ezNew);
	ezNew *= -opts_.dt;
	ezNew.Add(1.0, Ez_);

	// Update H.
	KxE_->Mult(ezNew, aux);
	MInv_->Mult(aux, hyNew);
	hyNew *= -opts_.dt;
	hyNew.Add(1.0, Hy_);

	Ez_ = ezNew;
	Hy_ = hyNew;
}

}