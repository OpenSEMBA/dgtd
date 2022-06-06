#ifndef _HEADER_MAXWELLDGTRACEINTEGRATOR
#define _HEADER_MAXWELLDGTRACEINTEGRATOR

#include "mfem.hpp"
#include "../../../general/forall.hpp"

namespace maxwell {

using namespace mfem;

class MaxwellWeakDGTraceIntegrator : public BilinearFormIntegrator
{

public:
	//When explicitly undeclared, rho = 1.0;
	MaxwellWeakDGTraceIntegrator(VectorCoefficient& u_, double a)
	{
		rho = NULL; u = &u_; alpha = a; beta = 0.5 * a;
	}

	MaxwellWeakDGTraceIntegrator(VectorCoefficient& u_, double a, double b)
	{
		rho = NULL; u = &u_; alpha = a; beta = b;
	}

	MaxwellWeakDGTraceIntegrator(Coefficient& rho_, VectorCoefficient& u_,
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

class MaxwellStrongDGTraceIntegrator : public TransposeIntegrator
{

public:
	MaxwellStrongDGTraceIntegrator(VectorCoefficient& u, double a)
		: TransposeIntegrator(new MaxwellWeakDGTraceIntegrator(u, -a, 0.5 * a)) { }

	MaxwellStrongDGTraceIntegrator(VectorCoefficient& u, double a, double b)
		: TransposeIntegrator(new MaxwellWeakDGTraceIntegrator(u, -a, b)) { }

	MaxwellStrongDGTraceIntegrator(Coefficient& rho, VectorCoefficient& u,
		double a, double b)
		: TransposeIntegrator(new MaxwellWeakDGTraceIntegrator(rho, u, -a, b)) { }

};
}

#endif
