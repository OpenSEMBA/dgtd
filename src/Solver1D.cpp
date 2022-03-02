#include "Solver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace Maxwell {

Solver::Solver(const Options& opts, const Mesh& mesh)
{

	checkOptionsAreValid(opts, mesh);

	mesh_ = mfem::Mesh(mesh, true);
	opts_ = opts;

	fec_ = std::make_unique<DG_FECollection>(opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
	fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	fecH1_ = std::make_unique<H1_FECollection>(opts_.order, mesh_.Dimension());
	fesH1_ = std::make_unique<FiniteElementSpace>(&mesh_, fecH1_.get());

	MInv_ = buildMassMatrix();
	Kx_ = buildDerivativeOperator(X);

	ez_.SetSpace(fes_.get());
	ez_.ProjectCoefficient(ConstantCoefficient(0.0));

	hx_.SetSpace(fes_.get());
	hx_.ProjectCoefficient(ConstantCoefficient(0.0));

	initializeParaviewData();

}
/*This test is a WIP.*/

int order = 2;
const int dimension = 1;
FiniteElementCollection* fec;
FiniteElementSpace* fes;

Mesh mesh = Mesh::MakeCartesian1D(100);
mesh.GetBoundingBox(
	AnalyticalFunctions::meshBoundingBoxMin,
	AnalyticalFunctions::meshBoundingBoxMax);
fec = new DG_FECollection(order, dimension, BasisType::GaussLobatto);
fes = new FiniteElementSpace(&mesh, fec);

BilinearForm massMatrixInv(fes);
massMatrixInv.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
massMatrixInv.Assemble();
massMatrixInv.Finalize();

ConstantCoefficient one(1.0);
BilinearForm K(fes);
K.AddDomainIntegrator(
	new TransposeIntegrator(
		new DerivativeIntegrator(one, 0)));

std::vector<VectorConstantCoefficient> n = {
	VectorConstantCoefficient(Vector({1.0})),
};

K.AddInteriorFaceIntegrator(
	new DGTraceIntegrator(n[0], -1.0, 0.0));
K.AddBdrFaceIntegrator(
	new DGTraceIntegrator(n[0], -1.0, 0.0));

K.Assemble();
K.Finalize();

GridFunction Ez(fes);
Ez.ProjectCoefficient(FunctionCoefficient(AnalyticalFunctions::gaussianFunction1D));
GridFunction Hy(fes);
Hy.ProjectCoefficient(ConstantCoefficient(0.0));

std::unique_ptr<ParaViewDataCollection> pd = NULL;
pd = std::make_unique<ParaViewDataCollection>("MaxwellView1D", &mesh);
pd->SetPrefixPath("ParaView");
pd->RegisterField("ez", &Ez);
pd->RegisterField("hy", &Hy);
pd->SetDataFormat(VTKFormat::BINARY);
pd->SetHighOrderOutput(true);

double time = 0.0;
bool done = false;

Vector aux(fes->GetVSize());
Vector ezNew(fes->GetVSize());
Vector hyNew(fes->GetVSize());

pd->SetCycle(0);
pd->SetTime(0.0);
pd->Save();

for (int cycle = 0; !done;)
{

	// Update E.
	K.Mult(Hy, aux);
	massMatrixInv.Mult(aux, ezNew);
	ezNew *= dt;
	ezNew.Add(1.0, Ez);

	// Update H.
	K.Mult(Ez, aux);
	massMatrixInv.Mult(aux, hyNew);
	hyNew *= dt;
	hyNew.Add(1.0, Hy);

	Ez = ezNew;
	Hy = hyNew;

	time += dt;
	cycle++;

	done = (time >= t_final - 1e-8 * dt);

	if (done || cycle % vis_steps == 0) {
		pd->SetCycle(cycle);
		pd->SetTime(time);
		pd->Save();
	}
}


/*This test is a WIP.*/

int order = 2;
const int dimension = 1;
FiniteElementCollection* fec;
FiniteElementSpace* fes;

Mesh mesh = Mesh::MakeCartesian1D(100);
mesh.GetBoundingBox(
	AnalyticalFunctions::meshBoundingBoxMin,
	AnalyticalFunctions::meshBoundingBoxMax);
fec = new DG_FECollection(order, dimension, BasisType::GaussLobatto);
fes = new FiniteElementSpace(&mesh, fec);

BilinearForm massMatrixInv(fes);
massMatrixInv.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
massMatrixInv.Assemble();
massMatrixInv.Finalize();

ConstantCoefficient one(1.0);
BilinearForm K(fes);
K.AddDomainIntegrator(
	new TransposeIntegrator(
		new DerivativeIntegrator(one, 0)));

std::vector<VectorConstantCoefficient> n = {
	VectorConstantCoefficient(Vector({1.0})),
};

K.AddInteriorFaceIntegrator(
	new DGTraceIntegrator(n[0], -1.0, 0.0));
K.AddBdrFaceIntegrator(
	new DGTraceIntegrator(n[0], -1.0, 0.0));

K.Assemble();
K.Finalize();

GridFunction Ez(fes);
Ez.ProjectCoefficient(FunctionCoefficient(AnalyticalFunctions::gaussianFunction1D));
GridFunction Hy(fes);
Hy.ProjectCoefficient(ConstantCoefficient(0.0));

std::unique_ptr<ParaViewDataCollection> pd = NULL;
pd = std::make_unique<ParaViewDataCollection>("MaxwellView1D", &mesh);
pd->SetPrefixPath("ParaView");
pd->RegisterField("ez", &Ez);
pd->RegisterField("hy", &Hy);
pd->SetDataFormat(VTKFormat::BINARY);
pd->SetHighOrderOutput(true);

double time = 0.0;
bool done = false;
double dt = 1e-4;
double t_final = 1000 * dt;
int vis_steps = 100;

Vector aux(fes->GetVSize());
Vector ezNew(fes->GetVSize());
Vector hyNew(fes->GetVSize());

pd->SetCycle(0);
pd->SetTime(0.0);
pd->Save();

for (int cycle = 0; !done;)
{

	// Update E.
	K.Mult(Hy, aux);
	massMatrixInv.Mult(aux, ezNew);
	ezNew *= dt;
	ezNew.Add(1.0, Ez);

	// Update H.
	K.Mult(Ez, aux);
	massMatrixInv.Mult(aux, hyNew);
	hyNew *= dt;
	hyNew.Add(1.0, Hy);

	Ez = ezNew;
	Hy = hyNew;

	time += dt;
	cycle++;

	done = (time >= t_final - 1e-8 * dt);

	if (done || cycle % vis_steps == 0) {
		pd->SetCycle(cycle);
		pd->SetTime(time);
		pd->Save();
	}
}


}
}