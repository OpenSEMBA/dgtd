#include "LinearIntegrators.h"
#include "IntegratorFunctions.h"

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

        double un{ vu * nor.operator*=(-1.0)};
        double w{ beta_ * un * ip.weight };
        elvect.Add(w, shape1_);
    }
}

void BoundaryDGJumpIntegrator::AssembleRHSElementVect(
    const FiniteElement& el1, const FiniteElement& el2, FaceElementTransformations& Tr, Vector& elvect)
{

    int dim{ el1.GetDim() };
    double vu1_data[3], nor1_data[3], vu2_data[3], nor2_data[3];
    Vector vu1(vu1_data, dim), nor1(nor1_data, dim), vu2(vu2_data, dim), nor2(nor2_data, dim);

    const IntegrationRule* ir = IntRule;
    if (ir == NULL)
    {
        int order;
        if (Tr.Elem2No >= 0) {
            order = (std::min(Tr.Elem1->OrderW(), Tr.Elem2->OrderW()) +
                2 * std::max(el1.GetOrder(), el2.GetOrder()));
        }
        else {
            order = Tr.Elem1->OrderW() + 2 * el1.GetOrder();
        }
        if (el1.Space() == FunctionSpace::Pk) {
            order++;
        }
        ir = &IntRules.Get(Tr.GetGeometryType(), order);

    }

    int ndof1 = el1.GetDof();
    int ndof2 = setNeighbourNDoF(el2, Tr);

    if (ndof1 > ndof2) {
        elvect.SetSize(ndof1);
    }
    else {
        elvect.SetSize(ndof2);
    }
    elvect = 0.0;

    for (int p = 0; p < ir->GetNPoints(); p++)
    {
        const IntegrationPoint& ip = ir->IntPoint(p);

        // Set the integration point in the face and the neighboring element
        Tr.SetAllIntPoints(&ip);

        // Access the neighboring element's integration point
        const IntegrationPoint& eip1 = Tr.GetElement1IntPoint();
        const IntegrationPoint& eip2 = Tr.GetElement2IntPoint();


        // Use Tr.Elem1 transformation for u so that it matches the coefficient
        // used with the ConvectionIntegrator and/or the DGTraceIntegrator.
        u_->Eval(vu1, *Tr.Elem1, eip1);
        u_->Eval(vu2, *Tr.Elem2, eip2);

        if (dim == 1)
        {
            nor1(0) = 2 * eip1.x - 1.0;
            nor2(0) = 2 * eip2.x - 1.0;
        }
        else
        {
            Vector ortho(dim);
            CalcOrtho(Tr.Jacobian(), ortho);
            nor1 = ortho.operator/=(Tr.Weight());
            nor2 = ortho.operator/=(Tr.Weight());
        }

        double un1{ vu1 * nor1.operator*=(-1.0) }, w1{ beta_ * un1 * ip.weight };
        double un2{ vu2 * nor2.operator*=(-1.0) }, w2{ beta_ * un2 * ip.weight };

        shape1_.SetSize(ndof1);
        shape2_.SetSize(ndof2);
        el1.CalcShape(eip1, shape1_);
        el2.CalcShape(eip2, shape2_);

        elvect.Add(w1, shape1_);
        elvect.Add(w2, shape2_);
    }
}

}
}