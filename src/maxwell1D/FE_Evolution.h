#pragma once

#include "mfem.hpp"

namespace Maxwell1D {

using namespace mfem;

class FE_Evolution : public TimeDependentOperator
{
public:
	FE_Evolution(
		BilinearForm& MInv,  
		BilinearForm& kTerm);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:
	BilinearForm &MInv_, &kTerm_;
	Vector z_;
};

FE_Evolution::FE_Evolution(BilinearForm& MInv, BilinearForm& kTerm): 
	TimeDependentOperator(MInv.Height()), 
	MInv_(MInv), 
	kTerm_(kTerm), 
	z_(MInv.Height())
{}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	Vector aux(MInv_.Height());

	// Update other term.
	kTerm_.Mult(x, aux);
	MInv_.Mult(aux, y);

	//// Update E.
	//Kh_.Mult(Hy, aux);
	//MInv_.Mult(aux, ezNew);
	// 
	//// Update H.
	//Ke_.Mult(Ez, aux);
	//MInv_.Mult(aux, hyNew);
}

}