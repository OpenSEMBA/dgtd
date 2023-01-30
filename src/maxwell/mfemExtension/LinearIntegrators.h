#pragma once

#include <mfem.hpp>

namespace maxwell {
namespace mfemExtension {

using namespace mfem;

/** Class for boundary integration of the linear form:
beta < (u.n) f, w >,
where f and u are given scalar and vector coefficients, respectively,
and w is the scalar test function. */
class BoundaryDGJumpIntegrator : public LinearFormIntegrator
{
private:
    Coefficient* f_;
    VectorCoefficient* u_;
    double beta_;

    Vector shape_;

public:
    BoundaryDGJumpIntegrator(Coefficient& f, VectorCoefficient& u,
        double beta)
    {
        f_ = &f; u_ = &u; beta_ = beta;
    }

    virtual void AssembleRHSElementVect(const FiniteElement& el,
        FaceElementTransformations& Tr,
        Vector& elvect);

    using LinearFormIntegrator::AssembleRHSElementVect;
};
}
}