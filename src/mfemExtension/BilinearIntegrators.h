#pragma once

#include <mfem.hpp>

namespace maxwell{
namespace mfemExtension {

using namespace mfem;

using Direction = int;

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
class MaxwellDGZeroNormalJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGZeroNormalJumpIntegrator(double b)
	{
		beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	double beta;

private:
	Vector shape1_, shape2_;
};

class MaxwellDGOneNormalJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGOneNormalJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
	{
		dir = dirTerms; beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	std::vector<Direction> dir;
	double beta;

private:
	Vector shape1_, shape2_;
};

class MaxwellDGTwoNormalJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGTwoNormalJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
	{
		dir = dirTerms; beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	std::vector<Direction> dir;
	double beta;

private:
	Vector shape1_, shape2_;
};

class MaxwellDGTraceJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGTraceJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
	{
		dir = dirTerms; beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
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

class MaxwellDGInteriorJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGInteriorJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
	{
		dir = dirTerms; beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
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

class MaxwellSMAJumpIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellSMAJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
	{
		dir = dirTerms; beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
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

class MaxwellDGFluxTotalFieldIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGFluxTotalFieldIntegrator(const std::vector<Direction>& dirTerms, double coeff, double b)
	{
		dir = dirTerms; TFSFCoeff_ = coeff; beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	std::vector<Direction> dir;
	double beta;
	int dim;

private:
	double TFSFCoeff_;
	Vector shape1_, shape2_;
};

class MaxwellDGPenaltyTotalFieldIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGPenaltyTotalFieldIntegrator(const std::vector<Direction>& dirTerms, const double coeff, double b)
	{
		dir = dirTerms; TFSFCoeff_ = coeff; beta = b;
	}

	void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	std::vector<Direction> dir;
	double beta;
	int dim;

private:
	double TFSFCoeff_;
	Vector shape1_, shape2_;
};

class TotalFieldScatteredFieldIntegrator : public BilinearFormIntegrator
{

public:
	TotalFieldScatteredFieldIntegrator(Coefficient& q, double b) :
		coeff(&q),  beta(b) {}

	void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

	void AssembleFaceMatrix(const FiniteElement& el1,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	Coefficient* coeff;
	double beta;
	int dim;

private:
	Vector shape1_;
};

}
}