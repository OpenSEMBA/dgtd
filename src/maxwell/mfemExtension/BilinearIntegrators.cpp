#include "BilinearIntegrators.h"

namespace maxwell {
namespace mfemExtension {

const IntegrationRule* setIntegrationRule(
    const FiniteElement& el1,
    const FiniteElement& el2,
    FaceElementTransformations& Trans)
{

    int order;
    // Assuming order(u)==order(mesh)
    if (Trans.Elem2No >= 0) {
        order = (std::min(Trans.Elem1->OrderW(), Trans.Elem2->OrderW()) +
            2 * std::max(el1.GetOrder(), el2.GetOrder()));
    }
    else {
        order = Trans.Elem1->OrderW() + 2 * el1.GetOrder();
    }
    if (el1.Space() == FunctionSpace::Pk) {
        order++;
    }
    return &IntRules.Get(Trans.GetGeometryType(), order);
}


const Vector setNormalVector1D(const int dim,
    const IntegrationPoint& eip1
)
{
    Vector res(dim);
    res(0) = 2 * eip1.x - 1.0;
    return res;
}

const Vector setNormalVector(const int dim,
    const IntegrationPoint& eip1,
    FaceElementTransformations& Trans
)
{
    Vector res(dim);
    CalcOrtho(Trans.Jacobian(), res);
    return res;
}

void buildFaceMatrix(double w, int ndofA, int ndofB, int offsetRow, int offsetCol,
    Vector shapeA, Vector shapeB, DenseMatrix& elmat) {
    for (int i = 0; i < ndofA; i++) {
        for (int j = 0; j < ndofB; j++)
        {
            double nonDiag = +1.0;
            if (offsetRow != offsetCol) {
                nonDiag = -1.0;
            }
            elmat(i + offsetRow, j + offsetCol) += nonDiag * w * shapeA(i) * shapeB(j);
        }
    }
}

const int setNeighbourNDoF(const FiniteElement& el2, FaceElementTransformations& Trans)
{
    if (Trans.Elem2No >= 0) {
       return el2.GetDof();
    }
    else {
       return 0;
    }
}

const double buildNormalTerm(const Vector& nor, const Direction& dir)
{
    std::vector<double> res{0.0, 0.0, 0.0};
    for (int i = 0; i < nor.Size(); i++) {
        res.at(i) += nor.Elem(i);
    }
    return res.at(dir);
}

const Vector calculateNormal(const FiniteElement& el, const IntegrationPoint& eip, FaceElementTransformations& Trans)
{
    Vector res(1);
    switch (el.GetDim()) {
    case 1:
        res.SetSize(1);
        res = setNormalVector1D(el.GetDim(), eip);
        break;
    default:
        res.SetSize(el.GetDim());
        res = setNormalVector(el.GetDim(), eip, Trans);
        break;
    }
    return res;
}


double calculateBetaTerm(Vector& nor, std::vector<Direction>& dir, const double beta)
{
    switch (dir.size()) {
    case 0:
        return beta; //[v] = (v1-v2)
    case 1:
    {
        double nIn = buildNormalTerm(nor, dir.at(0));
        return beta * nIn; //nIn * [v] = nIn * (v1-v2)
    }
    case 2:
    {
        double nIn = buildNormalTerm(nor, dir.at(0));
        double nOut = buildNormalTerm(nor, dir.at(1));
        return beta * nIn * nOut; //(nIn * [v]) * nOut = nIn * (v1-v2) * nOut
    }
    default:
        throw std::exception("Incorrect dimensions for dirTerms vector.");
    }
}

/*########################## MDG START ##########################*/
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

        a = 0.5 * alpha * un; // 0.5 * n * {v} = n * (v1+v2)/2
        b = beta * fabs(un); //|n| * [v] = |n| * (v1-v2)
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

void MaxwellDGTraceJumpIntegrator::AssembleFaceMatrix(const FiniteElement& el1,
    const FiniteElement& el2,
    FaceElementTransformations& Trans,
    DenseMatrix& elmat)
{   

    int ndof1 = el1.GetDof();
    int ndof2 = setNeighbourNDoF(el2, Trans);

    shape1_.SetSize(ndof1);
    shape2_.SetSize(ndof2);
    elmat.SetSize(ndof1 + ndof2);
    elmat = 0.0;

    const IntegrationRule* ir = IntRule;
    if (ir == NULL)
    {
        ir = setIntegrationRule(el1, el2, Trans);
    }

    for (int p = 0; p < ir->GetNPoints(); p++)
    {
        const IntegrationPoint& ip = ir->IntPoint(p);

        Trans.SetAllIntPoints(&ip);

        const IntegrationPoint& eip1 = Trans.GetElement1IntPoint();
        const IntegrationPoint& eip2 = Trans.GetElement2IntPoint();

        Vector nor = calculateNormal(el1, eip1, Trans);

        double b = calculateBetaTerm(nor, dir, beta);

        el1.CalcShape(eip1, shape1_);
        double w = ip.weight * b;
        if (w != 0.0) {
            buildFaceMatrix    (w, ndof1, ndof1,     0,     0, shape1_ , shape1_, elmat);
        }

        if (ndof2) {

            el2.CalcShape(eip2, shape2_);

            if (w != 0.0) {
                buildFaceMatrix(w, ndof1, ndof2, 0, ndof1, shape1_, shape2_, elmat);
                switch (dir.size()) {
                case 0:
                    buildFaceMatrix(w, ndof2, ndof1, ndof1,     0, shape2_, shape1_, elmat);
                    buildFaceMatrix(w, ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);
                    break;
                default:
                    buildFaceMatrix(-w, ndof2, ndof1, ndof1,     0, shape2_, shape1_, elmat);
                    buildFaceMatrix(-w, ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);
                    break;
                }
            }
        }
    }
}

void HesthavenDerivativeIntegrator::AssembleElementMatrix2(
    const FiniteElement& trial_fe,
    const FiniteElement& test_fe,
    ElementTransformation& Trans,
    DenseMatrix& elmat)
{
    int dim = trial_fe.GetDim();
    int trial_nd = trial_fe.GetDof();
    int test_nd = test_fe.GetDof();
    int spaceDim = Trans.GetSpaceDim();

    int i, l;
    double det;

    elmat.SetSize(trial_nd, test_nd);
    dshape.SetSize(test_nd, dim);
    dshapedxt.SetSize(test_nd, spaceDim);
    dshapedxi.SetSize(test_nd);
    invdfdx.SetSize(dim, spaceDim);
    shape.SetSize(trial_nd);

    const IntegrationRule* ir = IntRule;
    if (ir == NULL)
    {
        int order;
        if (trial_fe.Space() == FunctionSpace::Pk)
        {
            order = trial_fe.GetOrder() + test_fe.GetOrder() - 1;
        }
        else
        {
            order = trial_fe.GetOrder() + test_fe.GetOrder() + dim;
        }

        if (trial_fe.Space() == FunctionSpace::rQk)
        {
            ir = &RefinedIntRules.Get(trial_fe.GetGeomType(), order);
        }
        else
        {
            ir = &IntRules.Get(trial_fe.GetGeomType(), order);
        }
    }

    elmat = 0.0;
    for (i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint& ip = ir->IntPoint(i);

        
        test_fe.CalcDShape(ip, dshape);

        Trans.SetIntPoint(&ip);
        CalcInverse(Trans.Jacobian(), invdfdx);
        det = Trans.Weight();
        Mult(dshape, invdfdx, dshapedxt);
        for (l = 0; l < test_nd; l++)
        {
            dshapedxi(l) = dshapedxt(l, xi);
        }

        trial_fe.CalcShape(ip, shape);

        shape *= Q->Eval(Trans, ip) * det * ip.weight;
        AddMultVWt(dshapedxi, shape, elmat);
    }
}



}
}

