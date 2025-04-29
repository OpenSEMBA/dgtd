#include "LinearIntegrators.h"

namespace maxwell {
namespace mfemExtension {
    
using namespace mfem;

void BoundaryDGJumpIntegrator::AssembleRHSElementVect(
    const FiniteElement& el, ElementTransformation& Tr, Vector& elvect)
{
    mfem_error("BoundaryDGJumpIntegrator::AssembleRHSElementVect\n"
        "  is not implemented as boundary integrator!\n"
        "  Use LinearForm::AddBdrFaceIntegrator instead of\n"
        "  LinearForm::AddBoundaryIntegrator.");
}

void BoundaryDGJumpIntegrator::AssembleRHSElementVect(
    const FiniteElement& el, FaceElementTransformations& Tr, Vector& elvect)
{
        
    int dim{ el.GetDim() }; 
    double vu_data[3], nor_data[3];
    Vector vu(vu_data, dim), nor(nor_data, dim);

    const IntegrationRule* ir = IntRule;
    if (ir == NULL)
    {
        // Assuming order(u)==order(mesh)
        int order = Tr.Elem1->OrderW() + 2 * el.GetOrder();
        if (el.Space() == FunctionSpace::Pk)
        {
            order++;
        }
        ir = &IntRules.Get(Tr.GetGeometryType(), order);
    }

    int ndof{ el.GetDof() };
    shape1_.SetSize(ndof);
    elvect.SetSize(ndof);
    elvect = 0.0;

    for (int p = 0; p < ir->GetNPoints(); p++)
    {
        const IntegrationPoint& ip = ir->IntPoint(p);

        // Set the integration point in the face and the neighboring element
        Tr.SetAllIntPoints(&ip);

        // Access the neighboring element's integration point
        const IntegrationPoint& eip = Tr.GetElement1IntPoint();
        el.CalcShape(eip, shape1_);

        // Use Tr.Elem1 transformation for u so that it matches the coefficient
        // used with the ConvectionIntegrator and/or the DGTraceIntegrator.
        u_->Eval(vu, *Tr.Elem1, eip);

        if (dim == 1)
        {
            nor(0) = 2 * eip.x - 1.0;
        }
        else
        {
            Vector ortho(dim);
            CalcOrtho(Tr.Jacobian(), ortho);
            nor = ortho.operator/=(Tr.Weight());
        }

        double un{ vu * nor};
        double w{ beta_ * un * ip.weight };
        elvect.Add(w, shape1_);
    }
}

void BoundaryDGJumpIntegrator::AssembleRHSElementVect(
    const FiniteElement& el1, const FiniteElement& el2, FaceElementTransformations& Tr, Vector& elvect)
{

    double vu_data[3] , nor_data[3], 
           vu_data2[3], nor_data2[3];
    Vector vu (vu_data , el1.GetDim()), nor (nor_data , el1.GetDim()), 
           vu2(vu_data2, el2.GetDim()), nor2(nor_data2, el2.GetDim());

    const IntegrationRule* ir = IntRule;
    if (ir == NULL)
    {
        ir = setIntegrationRule(el1, el2, Tr);
    }

    elvect.SetSize(el1.GetDof() + el2.GetDof());
    elvect = 0.0;

    for (int p = 0; p < ir->GetNPoints(); p++)
    {
        const IntegrationPoint& ip = ir->IntPoint(p);

        // Set the integration point in the face and the neighboring element
        Tr.SetAllIntPoints(&ip);

        const IntegrationPoint& eip1 = Tr.GetElement1IntPoint();
        const IntegrationPoint& eip2 = Tr.GetElement2IntPoint();
        shape1_.SetSize(el1.GetDof());
        shape2_.SetSize(el2.GetDof());
        el1.CalcShape(eip1, shape1_);
        el2.CalcShape(eip2, shape2_);

        // Use Tr.Elem1 and Tr.Elem2 transformation for u so that it matches the coefficient used.
        u_->Eval(vu , *Tr.Elem1, eip1);
        u_->Eval(vu2, *Tr.Elem2, eip2);

        if (el1.GetDim() == 1) //Multiple dimensions are not possible in MFEM at the moment, any elX.GetDim() works.
        {
            nor(0)  = 2 * eip1.x - 1.0;
            nor2(0) = 2 * eip2.x - 1.0;
        }
        else
        {
            Vector ortho(el1.GetDim()); //Multiple dimensions are not possible in MFEM at the moment, any elX.GetDim() works.
            CalcOrtho(Tr.Jacobian(), ortho);
            nor = ortho.operator/=(Tr.Weight());
            nor2 = nor;
        }

        double un { vu * nor };
        double un2{ vu2 * nor2 };
        double w1{ beta_ * un  * ip.weight };
        double w2{ beta_ * un2 * ip.weight };
        
        for (int i = 0; i < shape1_.Size(); ++i) {
            elvect[i]                  = w1 * shape1_[i];
        }
        for (int i = 0; i < shape2_.Size(); ++i) {
            elvect[shape1_.Size() + i] = w2 * shape2_[i];
        }
    }
}

void FarFieldBdrFaceIntegrator::AssembleRHSElementVect(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, mfem::Vector& elvect)
{
    mfem_error("RCSBoundaryIntegrator::AssembleRHSElementVect\n"
        "  is not implemented as boundary integrator!\n"
        "  Use LinearForm::AddBdrFaceIntegrator instead of\n"
        "  LinearForm::AddBoundaryIntegrator.");
}

void FarFieldBdrFaceIntegrator::AssembleRHSElementVect(
    const mfem::FiniteElement& el1, const mfem::FiniteElement& el2, mfem::FaceElementTransformations& Tr, mfem::Vector& elvect)
{
    mfem_error("RCSBoundaryIntegrator::AssembleRHSElementVect\n"
        "  is not implemented for two element purposes!\n");
}

void FarFieldBdrFaceIntegrator::AssembleRHSElementVect(const FiniteElement& el, FaceElementTransformations& Tr, Vector& elvect) 
{
    // Initialise the shape and return vector, making the latter 0.0 as it will have things added onto it, 
    // and not overwritten.
    shape_.SetSize(el.GetDof());
    elvect.SetSize(el.GetDof());
    elvect = 0.0;

    // Construct or retrieve an integration rule for the appropriate reference element with the desired order of accuracy, 
    // taking into account the element's geometry type and an appropiate integration rule order.
    const IntegrationRule* ir = &IntRules.Get(Tr.GetGeometryType(), el.GetOrder() + Tr.Order());

    // Initialise vectors that will hold normal components.
    Vector inner_normal(3), normal(el.GetDim());

    // Loop over each quadrature point in the reference element
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        // Extract the current quadrature point from the integration rule
        const IntegrationPoint& ip = ir->IntPoint(i);

        // Prepare to evaluate the coordinate transformation at the current
        // quadrature point
        Tr.SetAllIntPoints(&ip);

        const IntegrationPoint& eip = Tr.GetElement1IntPoint();

        // We calculate the normal at the specified face, due to the problem 
        // we're solving and design choices, we invert said normal as we need it heading into the element.
        inner_normal = 0.0;
        normal = 0.0;
        CalcOrtho(Tr.Jacobian(), normal);

        for (auto i{ X }; i < normal.Size(); i++) {
            inner_normal[i] = -normal[i];
        }

        el.CalcShape(eip, shape_);

        // We evaluate the value of the coefficient at the specified point.
        auto coeff_eval{ c_.Eval(*Tr.Face, ip) };

        // Assemble the result of the calculation we wante to perform, that is
        // Weight of the IntegrationPoint * evaluation of the coefficient * normal on the specified direction / weight of the face surface, to make the normal vector unitary (Taflove p361 eq. 8.22a,b).
        auto val = ip.weight * coeff_eval * inner_normal[dir_];
        val /= Tr.Weight();

        elvect.Add(val, shape_);
    }

}
//
//void VectorNTFFBdrFaceIntegrator::AssembleRHSElementVect(
//    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, mfem::Vector& elvect)
//{
//    mfem_error("VectorNTFFBdrFaceIntegrator::AssembleRHSElementVect\n"
//        "  is not implemented as boundary integrator!\n"
//        "  Use LinearForm::AddBdrFaceIntegrator instead of\n"
//        "  LinearForm::AddBoundaryIntegrator.");
//}
//
//void VectorNTFFBdrFaceIntegrator::AssembleRHSElementVect(
//    const mfem::FiniteElement& el1, const mfem::FiniteElement& el2, mfem::FaceElementTransformations& Tr, mfem::Vector& elvect)
//{
//    mfem_error("VectorNTFFBdrFaceIntegrator::AssembleRHSElementVect\n"
//        "  is not implemented for two element purposes!\n");
//}
//
//const double calculateExponentialTerm(const double freq, const std::pair<double,double> phi_theta, const IntegrationPoint& ip, const bool isReal) 
//{
//    double res;
//    auto landa = physicalConstants::speedOfLight_SI / freq;
//    auto wavenumber = 2.0 * M_PI / landa;
//    auto rad_term = wavenumber * (ip.x * sin(phi_theta.second) * cos(phi_theta.first) + ip.y * sin(phi_theta.second) * sin(phi_theta.first) + ip.z * cos(phi_theta.second));
//    if (isReal) {
//        return std::cos(rad_term);
//    } else {
//        return std::sin(rad_term);
//    }
//}
//
//const Array<int> getNVDofsIDsForDim(const Direction& d, const double nvdofs)
//{
//    Array<int> res(nvdofs / 3);
//    for (auto i{ 0.0 }; i < nvdofs; i++) {
//        res[i] = i + d * nvdofs;
//    }
//    return res;
//}
//
//void VectorNTFFBdrFaceIntegrator::AssembleRHSElementVect(const FiniteElement& el, FaceElementTransformations& Tr, Vector& elvect)
//{
//
//    //const int vdim = 3;
//    //const int ndofs = el.GetDof();
//    //const int nvdofs = vdim * ndofs;
//
//    //// Initialise the shape and return vector. We will subvertly divide the problem into three components inside, then return the finished calculated result.
//    //// This means the shape vector will have ndofs size, but the return vector will be vdim times ndofs.
//    //shape_.SetSize(ndofs);
//    //elvect.SetSize(nvdofs);
//    //elvect = 0.0;
//
//    //// Construct or retrieve an integration rule for the appropriate reference element with the desired order of accuracy, 
//    //// taking into account the element's geometry type and an appropiate integration rule order.
//    //const IntegrationRule* ir = &IntRules.Get(Tr.GetGeometryType(), el.GetOrder() + Tr.Order());
//
//    //// Initialise vectors that will hold normal components.
//    //Vector inner_normal(3), normal(el.GetDim());
//    //Vector Ax, Ay, Az;
//    //A_.GetSubVector(getNVDofsIDsForDim(X, nvdofs), Ax);
//    //A_.GetSubVector(getNVDofsIDsForDim(Y, nvdofs), Ay);
//    //A_.GetSubVector(getNVDofsIDsForDim(Z, nvdofs), Az);
//
//    //// Loop over each quadrature point in the reference element
//    //for (int i = 0; i < ir->GetNPoints(); i++)
//    //{
//    //    // Extract the current quadrature point from the integration rule
//    //    const IntegrationPoint& ip = ir->IntPoint(i);
//
//    //    // We calculate the exponential related coefficient ahead
//
//    //    // Prepare to evaluate the coordinate transformation at the current
//    //    // quadrature point
//    //    Tr.SetAllIntPoints(&ip);
//
//    //    const IntegrationPoint& eip = Tr.GetElement1IntPoint();
//
//    //    // In the following lines we merely calculate the normal at the specified face, due to the problem 
//    //    // we're solving and design choices, we invert said normal as we need it heading into the element.
//    //    inner_normal = 0.0;
//    //    normal = 0.0;
//    //    CalcOrtho(Tr.Jacobian(), normal);
//
//    //    for (auto i{ X }; i < normal.Size(); i++) {
//    //        inner_normal[i] = -normal[i];
//    //    }
//
//    //    el.CalcShape(eip, shape_);
//
//    //    // We evaluate the value of the coefficient at the specified point.
//    //    auto coeff_eval{ calculateExponentialTerm(freq_, {phi_,theta_}, ip, isReal_) };
//
//    //    // Assemble the result of the calculation we wante to perform, that is
//    //    // Weight of the IntegrationPoint * evaluation of the coefficient * normal on the specified direction / weight of the face surface, to make the normal vector unitary (Taflove p361 eq. 8.22a,b).
//    //    auto val = ip.weight * coeff_eval * normal[dir_];
//    //    val /= Tr.Weight();
//
//    //    elvect.Add(val, shape_);
//    //}
//
//}

}
}