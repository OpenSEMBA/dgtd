#pragma once

#include <mfem.hpp>
#include <components/Types.h>

namespace maxwell{
namespace mfemExtension {

using namespace mfem;

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

class MaxwellDGDecoupledZeroNormalJumpIntegrator : public BilinearFormIntegrator
{
public:
    MaxwellDGDecoupledZeroNormalJumpIntegrator(double b)
    {
        beta = b;
    }

    void AssembleFaceMatrix(const FiniteElement& el1,
        const FiniteElement& el2,
        FaceElementTransformations& Trans,
        DenseMatrix& elmat) override;

protected:
    double beta;

private:
    Vector shape1_, shape2_;
};

class MaxwellDGDecoupledOneNormalJumpIntegrator : public BilinearFormIntegrator
{
public:
    MaxwellDGDecoupledOneNormalJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
    {
        dir = dirTerms; beta = b;
    }

    void AssembleFaceMatrix(const FiniteElement& el1,
        const FiniteElement& el2,
        FaceElementTransformations& Trans,
        DenseMatrix& elmat) override;

protected:
    std::vector<Direction> dir;
    double beta;

private:
    Vector shape1_, shape2_;
};

class MaxwellDGDecoupledTwoNormalJumpIntegrator : public BilinearFormIntegrator
{
public:
    MaxwellDGDecoupledTwoNormalJumpIntegrator(const std::vector<Direction>& dirTerms, double b)
    {
        dir = dirTerms; beta = b;
    }

    void AssembleFaceMatrix(const FiniteElement& el1,
        const FiniteElement& el2,
        FaceElementTransformations& Trans,
        DenseMatrix& elmat) override;

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
	MaxwellSMAJumpIntegrator(const std::vector<Direction>& dirTerms)
	{
		dir = dirTerms; beta = 1.0;
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

class MaxwellDGOneNormalTotalFieldIntegrator : public BilinearFormIntegrator
{

public:
	MaxwellDGOneNormalTotalFieldIntegrator(const std::vector<Direction>& dirTerms, double b)
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
	TotalFieldScatteredFieldIntegrator(double b) :
		beta(b) {}

	void AssembleFaceMatrix(const FiniteElement& el1,
		const FiniteElement& el2,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

	void AssembleFaceMatrix(const FiniteElement& el1,
		FaceElementTransformations& Trans,
		DenseMatrix& elmat);

protected:
	double beta;

private:
	Vector shape1_;
};

class HesthavenFluxIntegrator : public BilinearFormIntegrator
{

public:
	HesthavenFluxIntegrator(double b)
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

/** Like MFEM's DerivativeIntegrator but computes adj(J)*grad_ref(phi)
    directly via CalcAdjugate, avoiding the CalcInverse/det(J) cancellation
    that loses precision on curved elements with small det(J). */
class AdjugateDerivativeIntegrator : public BilinearFormIntegrator
{
public:
	AdjugateDerivativeIntegrator(Coefficient &q, int i) : Q(&q), xi(i) { }

	void AssembleElementMatrix(const FiniteElement &el,
		ElementTransformation &Trans,
		DenseMatrix &elmat) override
	{
		AssembleElementMatrix2(el, el, Trans, elmat);
	}

	void AssembleElementMatrix2(const FiniteElement &trial_fe,
		const FiniteElement &test_fe,
		ElementTransformation &Trans,
		DenseMatrix &elmat) override;

private:
	int xi;
	Coefficient *Q;
	DenseMatrix dshape_, adjJ_, dshape_adj_;
	Vector dshapedxi_, shape_;
};

}
}