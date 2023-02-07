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

}
}