#include <gtest/gtest.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <mfem.hpp>

#include "TestUtils.h"

using namespace mfem;

double xPositionFunction(const Vector& pos)
{
	return pos[0];
};

double negFunction(const Vector& pos)
{

	double normalizedPos = 0.0;
	double leftBoundary = 0.0, rightBoundary = 5.0;
	double length = rightBoundary - leftBoundary;
	normalizedPos = (pos[0] - leftBoundary) / length;
	
	return -2 * normalizedPos;
}

double linearDummyFunction(const Vector& v, double time)
{
	return v[0] + time;
}

double RotatedGaussianFunction(const Vector& pos, double time) {
	return 1.0 *
		exp(
			-pow(pos[0] - (3.5 * cos(-M_PI / 4.0)) + pos[1] - (3.5 * -sin(-M_PI / 4.0)), 2.0) /
			(2.0 * pow(0.2, 2.0))
		);
}


class GridFunctionTest : public ::testing::Test {
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
		const double length = 1.0,
		const decltype(BasisType::GaussLobatto)& bt = BasisType::GaussLobatto)
	{
		mesh_ = Mesh::MakeCartesian1D(elements, length);
		fec_ = std::make_unique<DG_FECollection>(order, 1, bt);
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

TEST_F(GridFunctionTest, ElementIDFinder1D)
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

TEST_F(GridFunctionTest, IntegrationPointFinder1D)
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

TEST_F(GridFunctionTest, GetValuesAtPoints1D)
{
	auto m{ Mesh::MakeCartesian1D(4) };
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &m, &fec };

	DenseMatrix pointMat{ 
		{ 
			{0.0}, // Point at left boundary.
			{0.1}, // Point within an element in the mesh.
			{0.5}, // Point in the boundary between two elements.
			{1.0}, // Point at right boundary.
			{1.1}  // Point out of the mesh.
		} 
	};
	pointMat.Transpose();

	Array<int> elArray;
	Array<IntegrationPoint> ipArray;
	m.FindPoints(pointMat, elArray, ipArray);

	GridFunction gridFunction{ &fes };
	FunctionCoefficient fc{ xPositionFunction };
	gridFunction.ProjectCoefficient(fc);

	Array<int> expectedID{
		{
			0, // Id of first element (at left boundary).
			0, // Id of first element.
			1, // Lower Id of the neighbouring elements.
			3, // Id of last element (at right boundary).
			-1 // Negative for failure to find.
		}
	};
	for (int i = 0; i < elArray.Size(); i++) {
		EXPECT_EQ(expectedID[i], elArray[i]);
	}
	Vector expectedValues({ 0.0, 0.1, 0.5, 1.0, 999.0 });
	for (int i = 0; i < elArray.Size(); i++) {
		if (elArray[i] < 0) {
			continue;
		}
		EXPECT_EQ(gridFunction.GetValue(elArray[i], ipArray[i]), expectedValues[i]);
	}

}

TEST_F(GridFunctionTest, SetValuesAtPoints1D)
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
	ConstantCoefficient oneCoeff{ 1.0 };
	gridBase.ProjectCoefficient(oneCoeff);

	GridFunction adder(fes_.get()), subber(fes_.get());
	adder = 0.0; 
	subber = 0.0;
	ConstantCoefficient halfCoeff{ 0.5 };
	adder.ProjectCoefficient(halfCoeff, positiveDoF);
	ConstantCoefficient minusHalfCoeff{ -0.5 };
	subber.ProjectCoefficient(minusHalfCoeff, negativeDoF);

	gridBase.operator+=(adder);
	gridBase.operator+=(subber);

	EXPECT_EQ(0.5, gridBase.Elem(1));
	EXPECT_EQ(1.5, gridBase.Elem(2));
	EXPECT_EQ(1.5, gridBase.Elem(7));
	EXPECT_EQ(0.5, gridBase.Elem(8));
}

TEST_F(GridFunctionTest, TimeDependentGridFunction)
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

TEST_F(GridFunctionTest, ProjectBdrFunction)
{
	const auto order{ 1 };
	
	auto mesh{ Mesh::MakeCartesian1D(5, 5.0) };
	mesh.AddBdrPoint(2, 2);

	H1_FECollection fec{ order, 1, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &mesh, &fec };

	GridFunction f{ &fes };
	f = 0.0;
	FunctionCoefficient function(linearDummyFunction);
	Array<int> attr{ 2 };
	f.ProjectBdrCoefficient(function, attr);

	EXPECT_EQ(0.0, f[0]);
	EXPECT_EQ(2.0, f[2]);
}

TEST_F(GridFunctionTest, ProjectFunctionOnMesh)
{
	int order{ 4 };
	
	auto mesh{ 
		Mesh::LoadFromFile(
			(mfemMeshes2DFolder() + "severalrotatedquadsnormalised.mesh").c_str(), 1, 0
		) 
	};
	DG_FECollection fec{ order, 2, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &mesh, &fec };

	GridFunction proj{ &fes };
	FunctionCoefficient RotatedGaussian(RotatedGaussianFunction);
	proj.ProjectCoefficient(RotatedGaussian);

	DenseMatrix pointMat{
		{
			{3.5 * cos(-M_PI / 4.0), 3.5 * -sin(-M_PI / 4.0) } // Center of Gaussian.
		}
	};
	pointMat.Transpose();

	Array<int> elArray;
	Array<IntegrationPoint> ipArray;
	mesh.FindPoints(pointMat, elArray, ipArray);

	Vector expectedValue({ 1.0 });
	EXPECT_NEAR(proj.GetValue(elArray[0], ipArray[0]), expectedValue[0], 1e-5);
}

TEST_F(GridFunctionTest, ProjectBetweenDifferentSpaces)
{
	// Choose any mesh in 2D, either MFEM-Cartesian or any other format.
	// auto mesh{ Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true) };
	auto mesh{ Mesh::LoadFromFile("testData/mfemMeshes/2D/square-disc.mesh") };

	// Let us assume we have a L2 space, 2-dimensional, with vdim equal to 1, meaning
	// our components are individually separated into a GridFunction representing X,
	// another representing Y, and a third one that would represent Z if it were the case. 
	auto dgfec{ DG_FECollection(3, 2) };
	auto dgfes{ FiniteElementSpace(&mesh, &dgfec) };

	// Additionally, we create a L2 space with vdim 2, which we will use later to create a 
	// GridFunction and store the component values in a single object.
	auto dgfesv2{ FiniteElementSpace(&mesh, &dgfec, 2) };

	// For the sake of having non-zero GridFunctions, we initialise them to 1.0 in all the
	// degrees of freedom of the problem. As the vdim 2 GridFunction will be filled with the
	// previous GridFunctions there is no need to initialise its values.
	// While we could just create the vdim 2 GridFunction, we are 'simulating' that our 
	// components are initially separated, thus this step is needed to be closer to our
	// Solver code.
	GridFunction dg_x(&dgfes);
	dg_x = 1.0;
	GridFunction dg_y(&dgfes);
	dg_y = 1.0;
	GridFunction dg_gf(&dgfesv2);

	// The vdim 2 GridFunction has its values ordered byNODES, this means first it stores the X
	// components, then the Y components, and if it were vdim 3, the Z components. So, we fill
	// the vector accordingly with our individual GridFunctions.
	dg_gf.SetVector(dg_x, 0);
	dg_gf.SetVector(dg_y, dg_x.Size());

	// We create a VectorGridFunctionCoefficient from the vdim 2 GridFunction.
	VectorGridFunctionCoefficient dg_vgfc(&dg_gf);

	// We create a Nedelec Collection and Finite Element Space with vdim 1, 
	// through all of this process, the order and dimension between all the
	// created FEC is the same. As Nedelec spaces are vectorial by nature, 
	// it doesn't seem we need to initialise the FES with vdim 2, which could feel
	// like a natural choice.
	auto ndfec{ ND_FECollection(3, 2) };
	auto ndfes{ FiniteElementSpace(&mesh, &ndfec) };

	// We create a Nedelec GridFunction using the Nedelec FES. As we will be 
	// projecting a solution onto the GridFunction, it is not needed to set its
	// values to 0.0 
	GridFunction nd_gf(&ndfes);

	// We project the L2 VGFC onto the Nedelec GridFunction we have previously created, 
	// it doesn't seem we need to specify where goes what, even if the debugger tells us
	// that technically the Nedelec GridFunction is only half the size of the VGFC vector, 
	// it just works. 
	// I personally still do not understand where the information for the Y-Component goes.
	// It might have to do with the DegreesOfFreedom(DoF) being Vector(?)DegreesOfFreedom (VDoF)
	// in this case, thus storing the information at a lower level we cannot trivially see in the debugger.
	nd_gf.ProjectCoefficient(dg_vgfc);

	// This is the typical paraview exporting code for the GridFunctions and problem.
	//ParaViewDataCollection* pd = NULL;
	//pd = new ParaViewDataCollection("L2toND", &mesh);
	//pd->SetPrefixPath("SpaceSwapping");
	//pd->RegisterField("Galerkin Solution X", &dg_x);
	//pd->RegisterField("Galerkin Solution Y", &dg_y);
	//pd->RegisterField("Nedelec Solution", &nd_gf);
	//pd->SetLevelsOfDetail(1);
	//pd->SetDataFormat(VTKFormat::BINARY);
	//pd->SetHighOrderOutput(true);
	//pd->SetCycle(0);
	//pd->SetTime(0.0);
	//pd->Save();
	
	Vector nd_x, nd_y;
	nd_gf.GetVectorFieldNodalValues(nd_x, 1);
	nd_gf.GetVectorFieldNodalValues(nd_y, 2);

	// Simple condition to automatise the test. Strange as it is, the component argument in 
	// GetVectorFieldNodalValues starts at 1 - X, 2 - Y, and not at 0. Only to be lowered by 1,
	// inside the method.
	double error{ 1e-6 };
	for (auto i{ 0 }; i < nd_x.Size(); ++i) {
		EXPECT_NEAR(1.0, nd_x[i], error);
		EXPECT_NEAR(1.0, nd_y[i], error);
	}


}