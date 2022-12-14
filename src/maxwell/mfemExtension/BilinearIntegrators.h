#pragma once

#include <mfem.hpp>

#include "../../../general/forall.hpp"

namespace maxwell{
namespace mfemExtension {

using namespace mfem;

using Direction = int;

class MaxwellDGTraceIntegrator : public mfem::BilinearFormIntegrator {

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
	double alpha, beta;
	// PA extension
	Vector pa_data;
	const DofToQuad* maps;             ///< Not owned
	const FaceGeometricFactors* geom;  ///< Not owned
	int dim, nf, nq, dofs1D, quad1D;

private:
	Vector shape1_, shape2_;
};

/** Integrator for a specialised application of the DG form:
	beta < k [v], [w] >,
	where v and w are the trial and test variables, respectively.
	[v] is the jump such that [v]=(v1-v2) for the  face between elements 1 and 2.
	For boundary elements, v2=0. Coefficient k is a std::vector<Direction> that can
	have either one, two or no entries based on the following form:

	Empty vector:
	beta <  [v], [w]>
	One entry vector:
	beta <  n_dir_0 [v], [w] >
	Two entry vector:
	beta < (n_dir_0 [v]) * n_dir_1, [w] >

	One use case for this integrator is to discretize the individual Maxwell
	equations with a DG formulation.
	*/
class MaxwellDGTraceJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGTraceJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
	{
		dir = dirTerms; beta = b;
	}

	virtual void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);
	
protected:
	std::vector<Direction> dir;
	double beta;
	int dim;

private:
	Vector shape1_, shape2_;
};
/** Integrator for a specialised application of the Hesthaven DG form:
	(Q * u, dvdx),
	where u and v are the trial and test variables, respectively.

	One use case for this integrator is to discretize the individual Maxwell
	equations with the Hesthaven DG formulation.
	*/
class HesthavenDerivativeIntegrator : public BilinearFormIntegrator
{
protected:
	Coefficient* Q;

private:
	int xi;
	DenseMatrix dshape, dshapedxt, invdfdx;
	Vector shape, dshapedxi;

public:
	HesthavenDerivativeIntegrator(Coefficient& q, int i) : Q(&q), xi(i) { }
	virtual void AssembleElementMatrix(const FiniteElement& el,
		ElementTransformation& Trans,
		DenseMatrix& elmat)
	{
		AssembleElementMatrix2(el, el, Trans, elmat);
	}
	virtual void AssembleElementMatrix2(const FiniteElement& trial_fe,
		const FiniteElement& test_fe,
		ElementTransformation& Trans,
		DenseMatrix& elmat);
};

}
}