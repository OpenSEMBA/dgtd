#include "BilinearIntegrators.h"
#include "IntegratorFunctions.h"

namespace maxwell {
namespace mfemExtension {

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
        if (ndof2) {
            el2.CalcShape(eip2, shape2_);
        }
        double w = ip.weight * b * 0.5;
        if (w != 0.0) {
            switch (dir.size()) {
            case 0:
                buildFaceMatrix     ( w, ndof1, ndof1,     0,     0, shape1_, shape1_, elmat);//TL
                if (ndof2) {
                    buildFaceMatrix ( w, ndof1, ndof2,     0, ndof1, shape1_, shape2_, elmat);//TR
                    buildFaceMatrix ( w, ndof2, ndof1, ndof1,     0, shape2_, shape1_, elmat);//BL
                    buildFaceMatrix ( w, ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);//BR
                }
                break;
            case 1:
                buildFaceMatrix     ( w, ndof1, ndof1,     0,     0, shape1_, shape1_, elmat);//TL
                if (ndof2) {           
                    buildFaceMatrix ( w, ndof1, ndof2,     0, ndof1, shape1_, shape2_, elmat);//TR
                    buildFaceMatrix (-w, ndof2, ndof1, ndof1,     0, shape2_, shape1_, elmat);//BL
                    buildFaceMatrix (-w, ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);//BR
                }
                break;
            case 2:
                w /= Trans.Weight();
                w /= Trans.Weight();
                buildFaceMatrix    ( w, ndof1, ndof1,     0,     0, shape1_, shape1_, elmat);//TL
                if (ndof2) {          
                    buildFaceMatrix( w, ndof1, ndof2,     0, ndof1, shape1_, shape2_, elmat);//TR
                    buildFaceMatrix( w, ndof2, ndof1, ndof1,     0, shape2_, shape1_, elmat);//BL
                    buildFaceMatrix( w, ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);//BR
                }
                break;
            default:
                throw std::exception("Wrong direction size.");
            }
        }
    }
}

void MaxwellDGFluxTotalFieldIntegrator::AssembleFaceMatrix(const FiniteElement& el1,
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
        elmat = 0.0;
        if (w != 0.0) {
            if (TFSFCoeff_ > 0) {
                buildFaceMatrix(  w, ndof1, ndof1, 0, 0, shape1_, shape1_, elmat);
            }
            else {
                buildFaceMatrix(0.0, ndof1, ndof1, 0, 0, shape1_, shape1_, elmat);
            }
        }

        if (ndof2) {

            el2.CalcShape(eip2, shape2_);

            if (w != 0.0) {
                if (TFSFCoeff_ > 0) {
                    buildFaceMatrix(0.0, ndof1, ndof2,     0, ndof1, shape1_, shape2_, elmat);
                    buildFaceMatrix(0.0, ndof2, ndof1, ndof1,    0, shape2_, shape1_, elmat);
                    buildFaceMatrix(  w, ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);
                }
                else {
                    buildFaceMatrix(  w, ndof1, ndof2,     0, ndof1, shape1_, shape2_, elmat);
                    buildFaceMatrix(  w, ndof2, ndof1, ndof1,     0, shape2_, shape1_, elmat);
                    buildFaceMatrix(0.0, ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);
                }
            }
        }
    }
}

void MaxwellDGPenaltyTotalFieldIntegrator::AssembleFaceMatrix(const FiniteElement& el1,
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
        double w = ip.weight * b * TFSFCoeff_;
        elmat = 0.0;

        if (ndof2) {

            el2.CalcShape(eip2, shape2_);

            if (w != 0.0) {
                buildFaceMatrix(w,   ndof1, ndof2,     0, ndof1, shape1_, shape2_, elmat);
                buildFaceMatrix(w,   ndof2, ndof2, ndof1, ndof1, shape2_, shape2_, elmat);
            }
        }
    }
}



}
}

