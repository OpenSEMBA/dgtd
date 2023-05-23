#pragma once

#include <mfem.hpp>

#include "components/Types.h"

namespace maxwell {
namespace mfemExtension {
    
using namespace mfem;
const IntegrationRule* setIntegrationRule(const FiniteElement& el1, const FiniteElement& el2, FaceElementTransformations& Trans);

const int setNeighbourNDoF(const FiniteElement& el2, FaceElementTransformations& Trans);

const Vector setNormalVector1D(const int dim, const IntegrationPoint& eip1);
const Vector setNormalVector(const int dim, const IntegrationPoint& eip1, FaceElementTransformations& Trans);
const double buildNormalTerm(const Vector& nor, const Direction& dir);
const Vector calculateNormal(const FiniteElement& el, const IntegrationPoint& eip, FaceElementTransformations& Trans);

double calculateBetaTerm(Vector& nor, std::vector<Direction>& dir, const double beta);

void buildFaceMatrix(double w, int ndofA, int ndofB, int offsetRow, int offsetCol,
    Vector shapeA, Vector shapeB, DenseMatrix& elmat);

const Vector calculateSMANormal(const FiniteElement& el, const IntegrationPoint& eip, FaceElementTransformations& Trans);
const Vector setNormalSMAVector(const int dim, const IntegrationPoint& eip1, FaceElementTransformations& Trans);
}
}