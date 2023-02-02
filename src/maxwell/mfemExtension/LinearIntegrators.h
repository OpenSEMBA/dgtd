#pragma once

#include <mfem.hpp>

namespace maxwell {
namespace mfemExtension {

/** Class for boundary integration of the linear form:
beta < (u.n) f, w >,
where f and u are given scalar and vector coefficients 
and w is the scalar test function. */
class BoundaryDGJumpIntegrator : public mfem::LinearFormIntegrator
{
public:
    BoundaryDGJumpIntegrator(mfem::Coefficient& f, mfem::VectorCoefficient& u,
        double beta)
    {
        f_ = &f; u_ = &u; beta_ = beta;
    }

    void AssembleRHSElementVect(const mfem::FiniteElement& el,
        mfem::ElementTransformation& Tr,
        mfem::Vector& elvect);
    void AssembleRHSElementVect(const mfem::FiniteElement& el,
        mfem::FaceElementTransformations& Tr,
        mfem::Vector& elvect);

private:
    double beta_;
    mfem::VectorCoefficient* u_;
    mfem::Coefficient* f_;


    mfem::Vector shape_;

};
}
}