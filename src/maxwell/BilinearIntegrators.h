#ifndef _HEADER_MAXWELLDGTRACEINTEGRATOR
#define _HEADER_MAXWELLDGTRACEINTEGRATOR

#include "mfem.hpp"
#include "../../../general/forall.hpp"
#include "Types.h"

namespace maxwell {

using namespace mfem;

class MaxwellDGTraceIntegrator : public BilinearFormIntegrator
{

public:
	//When explicitly undeclared, rho = 1.0;
	MaxwellDGTraceIntegrator(const DisForm form, VectorCoefficient& u_, double a)
	{
		form_ = form; rho = NULL; u = &u_; alpha = a; beta = 0.5 * a;
	}

	MaxwellDGTraceIntegrator(const DisForm form, VectorCoefficient& u_, double a, double b)
	{
		form_ = form; rho = NULL; u = &u_; alpha = a; beta = b;
	}

	MaxwellDGTraceIntegrator(const DisForm form, Coefficient& rho_, VectorCoefficient& u_,
		double a, double b)
	{
		form_ = form; rho = &rho_; u = &u_; alpha = a; beta = b;
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
	DisForm form_;
};

class MaxwellDGTraceJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGTraceJumpIntegrator(std::vector<Direction>& dirTerms, double b)
	{
		dir = dirTerms; alpha = 0; beta = b;
	}

	virtual void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);
	
protected:
	std::vector<Direction> dir;
	double alpha, beta;
	int dim;

private:
	Vector shape1_, shape2_;

	const IntegrationRule* setIntegrationRule(const FiniteElement& el1, const FiniteElement& el2, FaceElementTransformations& Trans);
	const Vector setNormalVector(const int dim, const IntegrationPoint& eip1, FaceElementTransformations& Trans);
	void buildFaceMatrix(double w, int ndofA, int ndofB, int desvI, int desvJ,
		Vector shapeA, Vector shapeB, DenseMatrix& elmat);
	const int setNeighbourNDoF(const FiniteElement& el2, FaceElementTransformations& Trans);
	const double buildNormalTerm(const Vector& innerNor, const Direction& outerDir);
};
}

#endif
