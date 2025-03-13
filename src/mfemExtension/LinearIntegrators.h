#pragma once

#include <mfem.hpp>
#include <components/Types.h>
#include <mfemExtension/IntegratorFunctions.h>
#include <math/PhysicalConstants.h>

namespace maxwell {
namespace mfemExtension {

/** Class for boundary integration of the linear form:
beta < (u.n) f, w >,
where f and u are given scalar and vector coefficients 
and w is the scalar test function. */
class BoundaryDGJumpIntegrator : public mfem::LinearFormIntegrator
{
public:
    BoundaryDGJumpIntegrator(mfem::VectorCoefficient& u, double beta)
    {
        u_ = &u; beta_ = beta;
    }

    void AssembleRHSElementVect(const mfem::FiniteElement& el,
        mfem::ElementTransformation& Tr,
        mfem::Vector& elvect);
    void AssembleRHSElementVect(const mfem::FiniteElement& el,
        mfem::FaceElementTransformations& Tr,
        mfem::Vector& elvect);
    void AssembleRHSElementVect(const mfem::FiniteElement& el1,
        const mfem::FiniteElement& el2,
        mfem::FaceElementTransformations& Tr,
        mfem::Vector& elvect);

private:
    double beta_;
    mfem::VectorCoefficient* u_;


    mfem::Vector shape1_;
    mfem::Vector shape2_;

};

class FarFieldBdrFaceIntegrator : public mfem::LinearFormIntegrator
{
public:
    FarFieldBdrFaceIntegrator(mfem::Coefficient& c, const Direction& outputDir) :
        c_(c), dir_(outputDir) {}

    void AssembleRHSElementVect(const mfem::FiniteElement& el,
        mfem::ElementTransformation& Tr,
        mfem::Vector& elvect);

    void AssembleRHSElementVect(const mfem::FiniteElement& el1,
        const mfem::FiniteElement& el2,
        mfem::FaceElementTransformations& Tr,
        mfem::Vector& elvect);

    void AssembleRHSElementVect(const mfem::FiniteElement& el, 
        mfem::FaceElementTransformations& Tr, 
        mfem::Vector& elvect);

private:
    mfem::Coefficient& c_;
    Direction dir_;

    mfem::Vector shape_;

};
//
//class VectorNTFFBdrFaceIntegrator : public mfem::LinearFormIntegrator {
//public:
//    VectorNTFFBdrFaceIntegrator(const std::vector<std::complex<double>>& FE, const Vector& H, const double freq, const double phi, const double theta = 0.0, const bool isReal = true) :
//        A_(A), freq_(freq), phi_(phi), theta_(theta), isReal_(isReal) {}
//
//    void AssembleRHSElementVect(const mfem::FiniteElement& el,
//        mfem::ElementTransformation& Tr,
//        mfem::Vector& elvect);
//
//    void AssembleRHSElementVect(const mfem::FiniteElement& el1,
//        const mfem::FiniteElement& el2,
//        mfem::FaceElementTransformations& Tr,
//        mfem::Vector& elvect);
//
//    void AssembleRHSElementVect(const mfem::FiniteElement& el,
//        mfem::FaceElementTransformations& Tr,
//        mfem::Vector& elvect);
//
//private:
//    const Vector& A_;
//    double freq_;
//    double phi_;
//    double theta_;
//    bool isReal_;
//
//    mfem::Vector shape_;
//
//};
}
}