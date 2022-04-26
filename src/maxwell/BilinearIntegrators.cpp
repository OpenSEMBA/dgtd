#include "BilinearIntegrators.h"

namespace maxwell {

/*########################## MDG START ##########################*/
//Has alpha (Done?)
void MaxwellDGTraceIntegrator::AssembleFaceMatrix(const FiniteElement& el1,
    const FiniteElement& el2,
    FaceElementTransformations& Trans,
    DenseMatrix& elmat)
{
    int dim, ndof1, ndof2;

    double un, a, b, w;

    dim = el1.GetDim();
    ndof1 = el1.GetDof();
    Vector vu(dim), nor(dim);

    if (Trans.Elem2No >= 0)
    {
        ndof2 = el2.GetDof();
    }
    else
    {
        ndof2 = 0;
    }

    shape1_.SetSize(ndof1);
    shape2_.SetSize(ndof2);
    elmat.SetSize(ndof1 + ndof2);
    elmat = 0.0;

    const IntegrationRule* ir = IntRule;
    if (ir == NULL)
    {
        int order;
        // Assuming order(u)==order(mesh)
        if (Trans.Elem2No >= 0)
            order = (std::min(Trans.Elem1->OrderW(), Trans.Elem2->OrderW()) +
                2 * std::max(el1.GetOrder(), el2.GetOrder()));
        else
        {
            order = Trans.Elem1->OrderW() + 2 * el1.GetOrder();
        }
        if (el1.Space() == FunctionSpace::Pk)
        {
            order++;
        }
        ir = &IntRules.Get(Trans.GetGeometryType(), order);
    }

    for (int p = 0; p < ir->GetNPoints(); p++)
    {
        const IntegrationPoint& ip = ir->IntPoint(p);

        // Set the integration point in the face and the neighboring elements
        Trans.SetAllIntPoints(&ip);

        // Access the neighboring elements' integration points
        // Note: eip2 will only contain valid data if Elem2 exists
        const IntegrationPoint& eip1 = Trans.GetElement1IntPoint();
        const IntegrationPoint& eip2 = Trans.GetElement2IntPoint();

        el1.CalcShape(eip1, shape1_);

        u->Eval(vu, *Trans.Elem1, eip1);

        if (dim == 1)
        {
            nor(0) = 2 * eip1.x - 1.0;
        }
        else
        {
            CalcOrtho(Trans.Jacobian(), nor);
        }
        
        un = vu * nor;
        a = 0.5 * alpha * un;
        b = beta;
        // note: if |alpha/2|==|beta| then |a|==|b|, i.e. (a==b) or (a==-b)
        //       and therefore two blocks in the element matrix contribution
        //       (from the current quadrature point) are 0

        if (rho)
        {
            double rho_p;
            if (un >= 0.0 && ndof2)
            {
                rho_p = rho->Eval(*Trans.Elem2, eip2);
            }
            else
            {
                rho_p = rho->Eval(*Trans.Elem1, eip1);
            }
            a *= rho_p;
            b *= rho_p;
        }

        w = ip.weight * (a + b);
        if (w != 0.0)
        {
            for (int i = 0; i < ndof1; i++)
                for (int j = 0; j < ndof1; j++)
                {
                    elmat(i, j) += w * shape1_(i) * shape1_(j);
                }
        }

        if (ndof2)
        {
            el2.CalcShape(eip2, shape2_);

            if (w != 0.0)
                for (int i = 0; i < ndof2; i++)
                    for (int j = 0; j < ndof1; j++)
                    {
                        elmat(ndof1 + i, j) -= w * shape2_(i) * shape1_(j);
                    }

            w = ip.weight * (b - a);
            if (w != 0.0)
            {
                for (int i = 0; i < ndof2; i++)
                    for (int j = 0; j < ndof2; j++)
                    {
                        elmat(ndof1 + i, ndof1 + j) += w * shape2_(i) * shape2_(j);
                    }

                for (int i = 0; i < ndof1; i++)
                    for (int j = 0; j < ndof2; j++)
                    {
                        elmat(i, ndof1 + j) -= w * shape1_(i) * shape2_(j);
                    }
            }
        }
    }
}
}


