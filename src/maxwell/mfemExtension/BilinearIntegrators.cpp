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
    if (Trans.Elem2No >= 0) {
        Vector res(dim);
        CalcOrtho(Trans.Jacobian(), res);
        return res;
    }
    else {
        Vector ortho(dim);
        CalcOrtho(Trans.Jacobian(), ortho);
        return ortho.operator/=(Trans.Weight());
    }
}

void buildFaceMatrix(double w, int ndofA, int ndofB, int offsetRow, int offsetCol,
    Vector shapeA, Vector shapeB, DenseMatrix& elmat) {
    for (int i = 0; i < ndofA; i++) {
        for (int j = 0; j < ndofB; j++)
        {
            if (offsetRow == offsetCol) {
                elmat(i + offsetRow, j + offsetCol) += w * shapeA(i) * shapeB(j);
            }
            else {
                elmat(i + offsetRow, j + offsetCol) -= w * shapeA(i) * shapeB(j);
            }
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

}
}

