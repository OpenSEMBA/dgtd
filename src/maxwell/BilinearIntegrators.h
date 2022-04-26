#ifndef _HEADER_MAXWELLDGTRACEINTEGRATOR
#define _HEADER_MAXWELLDGTRACEINTEGRATOR

#include "mfem.hpp"
#include "../../../general/forall.hpp"

namespace maxwell {

using namespace mfem;

class MaxwellDGTraceIntegrator : public BilinearFormIntegrator
{

public:
	//When explicitly undeclared, rho = 1.0;
	MaxwellDGTraceIntegrator(VectorCoefficient& u_, double a)
	{
		rho = NULL; u = &u_; alpha = a; beta = 0.5 * a;
	}

	MaxwellDGTraceIntegrator(VectorCoefficient& u_, double a, double b)
	{
		rho = NULL; u = &u_; alpha = a; beta = b;
	}

	MaxwellDGTraceIntegrator(Coefficient& rho_, VectorCoefficient& u_,
		double a, double b)
	{
		rho = &rho_; u = &u_; alpha = a; beta = b;
	}

	virtual void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	Coefficient* rho;
	VectorCoefficient* u;
	double alpha, beta, gamma;
	// PA extension
	Vector pa_data;
	const DofToQuad* maps;             ///< Not owned
	const FaceGeometricFactors* geom;  ///< Not owned
	int dim, nf, nq, dofs1D, quad1D;

private:
	Vector shape1_, shape2_;
};
}


#endif
