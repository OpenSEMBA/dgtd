#include "gtest/gtest.h"

#include "GlobalFunctions.h"
#include <vector>
#include <mfem.hpp>

using namespace mfem;

double posFunction(const Vector& pos)
{
	double normalizedPos = 0.0;
	double leftBoundary = 0.0, rightBoundary = 5.0;
	double length = rightBoundary - leftBoundary;
	normalizedPos = (pos[0] - leftBoundary) / length;

	return 2 * normalizedPos;
};

double negFunction(const Vector& pos)
{

	double normalizedPos = 0.0;
	double leftBoundary = 0.0, rightBoundary = 5.0;
	double length = rightBoundary - leftBoundary;
	normalizedPos = (pos[0] - leftBoundary) / length;
	
	return -2 * normalizedPos;
}

class TestGridFunction : public ::testing::Test {
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

TEST_F(TestGridFunction, ElementIDFinder1D)
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

TEST_F(TestGridFunction, IntegrationPointFinder1D)
{
	setFES1D(2, 5, 5.0);

	DenseMatrix pointMat({ {0.3},{2.5},{4.0} });
	pointMat.Transpose();
	mesh_.FindPoints(pointMat, elArray_, ipArray_);

	IntegrationPoint a, b, c;
	a.Set(0.3, 0.0, 0.0, 1.0);
	b.Set(0.5, 0.0, 0.0, 1.0);
	c.Set(1.0, 0.0, 0.0, 1.0);
	Array<IntegrationPoint> expectedIP({ a,b,c });

	for (int i = 0; i < elArray_.Size(); i++) {
		EXPECT_EQ(expectedIP[i].x, ipArray_[i].x);
	}

}

TEST_F(TestGridFunction, DoFFinder1D)
{
	/*Create a 1D FES mesh with the following form

	   0      1      2      3      4
	|-----||-----||-----||-----||-----|
	0     12     34     56     78     9 
	
	Where the goal is to find DoF tags for nodes 2 and 7, 
	and 1 and 8 through Physical Position values, 
	and store them in Array<int> containers. */

	// OPENED ISSUE DUE TO FINDPOINTS POTENTIAL ISSUE

	setFES1D(1, 5, 5.0);

	/* Node positions for... 1        2        7        8 */
	DenseMatrix pointMat({ {1.0}, {1.00001}, {4.0}, {4.00001} });
	pointMat.Transpose();
	mesh_.FindPoints(pointMat, elArray_, ipArray_);

	GridFunction positions(fes_.get());
	mesh_.GetNodes(positions);










}

TEST_F(TestGridFunction, GetGFValuesAtPoints1D)
{
	setFES1D(2, 5, 5.0);
	
	DenseMatrix pointMat_({ {0.3},{2.5},{4.0} });
	pointMat_.Transpose();
	mesh_.FindPoints(pointMat_, elArray_, ipArray_);

	GridFunction gridFunction(fes_.get());
	gridFunction.ProjectCoefficient(FunctionCoefficient(posFunction));

	Vector expectedValues({ 0.12, 1.0, 1.6 });
	for (int i = 0; i < elArray_.Size(); i++) {
		EXPECT_EQ(gridFunction.GetValue(elArray_[i], ipArray_[i]), expectedValues[i]);
	}

}

TEST_F(TestGridFunction, SetGFValuesAtPoints1D)
{
	/*Create a 1D FES playground with the following form
	
	   0      1      2      3      4
	|-----||-----||-----||-----||-----|
	0     12     34     56     78     9
	      -+                   +-  
	
	Where the goal is to project certain Coefficients on a GridFunction
	which has assigned the FES we've built. The DoF with a + sign
	will have a positive addition to their initially projected
	value (0.5) and the - sign will have a negative
	value added in the form of (-0.5).
	
	This is a prelude to PlaneWave implementation.*/
	
	setFES1D(1, 5, 5.0);

	Array<int> positiveDoF{ {2,7} }, negativeDoF{ {1,8} };

	GridFunction gridBase(fes_.get());
	gridBase.ProjectCoefficient(ConstantCoefficient{ 1.0 });

	GridFunction adder(fes_.get()), subber(fes_.get());
	adder = 0.0; 
	subber = 0.0;
	adder.ProjectCoefficient(ConstantCoefficient{ 0.5 }, positiveDoF);
	subber.ProjectCoefficient(ConstantCoefficient{ -0.5 }, negativeDoF);

	gridBase.operator+=(adder);
	gridBase.operator+=(subber);

	EXPECT_EQ(0.5, gridBase.Elem(1));
	EXPECT_EQ(1.5, gridBase.Elem(2));
	EXPECT_EQ(1.5, gridBase.Elem(7));
	EXPECT_EQ(0.5, gridBase.Elem(8));
}

double linearDummyFunction(const Vector& v, double time)
{
	return v[0] + time;
}

TEST_F(TestGridFunction, TimeDependentGridFunction)
{
	setFES1D(1, 5, 5.0);
	GridFunction f{ fes_.get() };
	FunctionCoefficient function(linearDummyFunction);
	
	f.ProjectCoefficient(function);
	EXPECT_EQ(0.0, f[0]);

	for (auto time{ 0.0 }; time <= 5.0; time += 1.0) {
		function.SetTime(time);
		f.ProjectCoefficient(function);
	}
	EXPECT_EQ(5.0, f[0]);
}
	
TEST_F(TestGridFunction, ProjectBdrFunction)
{
	const auto order{ 1 };
	
	auto mesh{ Mesh::MakeCartesian1D(5, 5.0) };
	mesh.AddBdrPoint(2, 2);

	H1_FECollection fec{ order, 1, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &mesh, &fec };

	GridFunction f{ &fes };
	f = 0.0;
	FunctionCoefficient function(linearDummyFunction);
	f.ProjectBdrCoefficient(function, Array<int>{ 2 });

	EXPECT_EQ(0.0, f[0]);
	EXPECT_EQ(2.0, f[2]);
}

