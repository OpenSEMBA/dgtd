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

        // Access the neighboring element's integration point
        const IntegrationPoint& eip1 = Tr.GetElement1IntPoint();
        const IntegrationPoint& eip2 = Tr.GetElement2IntPoint();
        shape1_.SetSize(el1.GetDof());
        shape2_.SetSize(el2.GetDof());
        el1.CalcShape(eip1, shape1_);
        el2.CalcShape(eip2, shape2_);

        // Use Tr.Elem1 transformation for u so that it matches the coefficient used.
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

}
}