#include "gtest/gtest.h"

#include "GlobalFunctions.h"
#include <vector>
#include <mfem.hpp>

using namespace mfem;

double simpleLinearFunction(const Vector& pos)
{
	double normalizedPos = 0.0;
	double leftBoundary = 0.0, rightBoundary = 5.0;
	double length = rightBoundary - leftBoundary;
	normalizedPos = (pos[0] - leftBoundary) / length;

	return 2 * normalizedPos;
};

class TestPointFinder : public ::testing::Test {
protected:

	typedef std::size_t Direction;

	void SetUp() override
	{
		mesh_ = Mesh::MakeCartesian1D(1);
		fec_ = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void setFES1D(
		const int order,
		const int elements = 1,
		const double length = 1.0)
	{
		mesh_ = Mesh::MakeCartesian1D(elements, length);
		fec_ = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	}

	void setFES2D(
		const int order,
		const int xElem = 1,
		const int yElem = 1)
	{
		mesh_ = Mesh::MakeCartesian2D(xElem, yElem, Element::Type::QUADRILATERAL);
		fec_ = std::make_unique<DG_FECollection>(order, mesh_.Dimension(), BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	Array<int> elArray_;
	Array<IntegrationPoint> ipArray_;

	double tol_ = 1e-5;

};

TEST_F(TestPointFinder, ElementIDFinder1D)
{
	setFES1D(2, 5, 5.0);
	
	DenseMatrix pointMat_({{0.3},{2.5},{4.0}});
	pointMat_.Transpose();
	mesh_.FindPoints(pointMat_,elArray_,ipArray_);

	Array<int> expectedID({ 0, 2, 3 });

	for (int i = 0; i < elArray_.Size(); i++) {
		EXPECT_EQ(expectedID[i], elArray_[i]);
	}

}

TEST_F(TestPointFinder, IntegrationPointFinder1D)
{
	setFES1D(2, 5, 5.0);

	DenseMatrix pointMat_({ {0.3},{2.5},{4.0} });
	pointMat_.Transpose();
	mesh_.FindPoints(pointMat_, elArray_, ipArray_);

	IntegrationPoint a, b, c;
	a.Set(0.3, 0.0, 0.0, 1.0);
	b.Set(0.5, 0.0, 0.0, 1.0);
	c.Set(1.0, 0.0, 0.0, 1.0);
	Array<IntegrationPoint> expectedIP({ a,b,c });

	for (int i = 0; i < elArray_.Size(); i++) {
		EXPECT_EQ(expectedIP[i].x, ipArray_[i].x);
	}

}


TEST_F(TestPointFinder, GetGFValuesAtPoints1D)
{
	setFES1D(2, 5, 5.0);
	
	DenseMatrix pointMat_({ {0.3},{2.5},{4.0} });
	pointMat_.Transpose();
	mesh_.FindPoints(pointMat_, elArray_, ipArray_);

	GridFunction gridFunction(fes_.get());
	gridFunction.ProjectCoefficient(FunctionCoefficient(simpleLinearFunction));

	Vector expectedValues({ 0.12, 1.0, 1.6 });
	for (int i = 0; i < elArray_.Size(); i++) {
		EXPECT_EQ(gridFunction.GetValue(elArray_[i], ipArray_[i]), expectedValues[i]);
	}

}

TEST_F(TestPointFinder, SetGFValuesAtPoints1D)
{
	setFES1D(2, 5, 5.0);

	DenseMatrix pointMat_({ {0.3},{2.5},{4.0} });
	pointMat_.Transpose();
	mesh_.FindPoints(pointMat_, elArray_, ipArray_);

	GridFunction gridBase(fes_.get());
	gridBase.ProjectCoefficient(FunctionCoefficient(simpleLinearFunction));
	GridFunction gridAdder(fes_.get());
	gridAdder = 0.0;
	

}

