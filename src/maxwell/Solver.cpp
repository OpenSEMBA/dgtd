#include <fstream>
#include <iostream>
#include <algorithm>

#include "Solver.h"

using namespace mfem;

namespace maxwell {

Solver::Solver(const Model& model, const Probes& probes,
const Sources& sources, const Options& options) :

model_(model),
probes_(probes),
sources_(sources),
opts_(options),
mesh_(model_.getMesh())

{
fec_ = std::make_unique<DG_FECollection>(
opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

odeSolver_ = std::make_unique<RK4Solver>();

maxwellEvol_ = std::make_unique<FiniteElementEvolutionNoCond>(fes_.get(), opts_.evolutionOperatorOptions, model_);

sol_ = Vector(FiniteElementEvolutionNoCond::numberOfFieldComponents * FiniteElementEvolutionNoCond::numberOfMaxDimensions * fes_->GetNDofs());
sol_ = 0.0;

for (int d = X; d <= Z; d++) {
	E_[d].SetSpace(fes_.get());
	E_[d].SetData(sol_.GetData() + d*fes_->GetNDofs());
	H_[d].SetSpace(fes_.get());
	H_[d].SetData(sol_.GetData() + (d+3)*fes_->GetNDofs());
}

//if (source_) {
	setInitialField();
//}

if (probes_.paraview) {
	initializeParaviewData();
}
if (probes_.glvis) {
	//initializeGLVISData(); //TODO
}
if (probes_.extractDataAtPoints) {
	auto aux = Solver::buildElemAndIntegrationPointArrays(probes_.integPointMat);
	elemIds_ = aux.first;
	integPointSet_ = Solver::buildIntegrationPointsSet(aux.second);
	fieldToExtract_ = probes_.fieldToExtract;
}
}

void Solver::checkOptionsAreValid(const Options& opts)
{
	if ((opts.order < 0) ||
		(opts.t_final < 0) ||
		(opts.dt < 0)) {
		throw std::exception("Incorrect parameters in Options");
	}
}

void Solver::setInitialField()
{
	for (int i = 0; i < sources_.getSourceVector().size(); i++) {

		auto source = sources_.getSourceVector().at(i);
		std::function<double(const Position&)> f = std::bind(&Source::evalGaussianFunction1D, &source, std::placeholders::_1);

		switch (source.getFieldType()) {
		case FieldType::E:
			E_[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
			return;
		case FieldType::H:
			H_[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
			return;
		}
	}
}

const GridFunction& Solver::getFieldInDirection(const FieldType& ft, const Direction& d) const
{
	switch (ft) {
	case FieldType::E:
		return E_[d];
	case FieldType::H:
		return H_[d];
	}
}

const Vector& Solver::getMaterialProperties(const Material& mat) const
{
	return Vector({mat.getPermittivity(), mat.getPermeability(), mat.getImpedance(), mat.getConductance()});
}

std::pair<Array<int>, Array<IntegrationPoint>>& Solver::buildElemAndIntegrationPointArrays(DenseMatrix& physPoints)
{
	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	fes_->GetMesh()->FindPoints(physPoints, elemIdArray, integPointArray);
	return { std::make_pair(elemIdArray, integPointArray) };
}

const std::vector<std::vector<IntegrationPoint>>& Solver::buildIntegrationPointsSet(const Array<IntegrationPoint>& ipArray) const
{
	IntegrationPointsSet res; 
	for (int i = 0; i < ipArray.Size(); i++) {
		switch (fes_->GetMesh()->Dimension()) {
		case 1:
			res[i][X].Set1w(ipArray[i].x, 0.0);
			break;
		case 2:
			res[i][X].Set2(ipArray[i].x, 0.0);
			res[i][Y].Set2(0.0, ipArray[i].y);
			break;
		case 3:
			res[i][X].Set3(ipArray[i].x, 0.0, 0.0);
			res[i][Y].Set3(0.0, ipArray[i].y, 0.0);
			res[i][Z].Set3(0.0, 0.0, ipArray[i].z);
			break;
		}
	}
	return res;
}
const std::array<std::array<double, 3>, 3>& Solver::saveFieldAtPoints(const FieldType& ft)
{
	auto maxDim = fes_->GetMesh()->Dimension();
	std::array<std::array<double, 3>, 3> res{0.0};
	for (int i = 0; i < elemIds_.Size(); i++) {
		for (int dir = Direction::X; dir != maxDim; dir++) {
			Direction d = static_cast<Direction>(dir);
			switch (ft) {
			case FieldType::E:
				res[i][d] = E_[d].GetValue(elemIds_[i],integPointSet_[i][d]);
			case FieldType::H:
				res[i][d] = H_[d].GetValue(elemIds_[i],integPointSet_[i][d]);
			}
		}
	}
	return res;
}


void Solver::initializeParaviewData()
{
	pd_ = NULL;
	pd_ = std::make_unique<ParaViewDataCollection>("MaxwellView", &mesh_);
	pd_->SetPrefixPath("ParaView");
	pd_->RegisterField("Ex", &E_[X]);
	pd_->RegisterField("Ey", &E_[Y]);
	pd_->RegisterField("Ez", &E_[Z]);
	pd_->RegisterField("Hx", &H_[X]);
	pd_->RegisterField("Hy", &H_[Y]);
	pd_->RegisterField("Hz", &H_[Z]);
	pd_->SetLevelsOfDetail(opts_.order);
	pd_->SetDataFormat(VTKFormat::BINARY);
	opts_.order > 0 ? pd_->SetHighOrderOutput(true) : pd_->SetHighOrderOutput(false);
}

//void Solver1D::initializeGLVISData() //TODO
//{
//	char vishost[] = "localhost";
//	int  visport = 19916;
//	sout_.open(vishost, visport);
//	sout_.precision(probes_.precision);
//	sout_ << "solution\n" << mesh_ << E_;
//	sout_ << "pause\n";
//	sout_ << std::flush;
//	std::cout << "GLVis visualization paused."
//		<< " Press space (in the GLVis window) to resume it.\n";
//}

void Solver::storeInitialVisualizationValues()
{
	if (probes_.paraview) {
		pd_->SetCycle(0);
		pd_->SetTime(0.0);
		pd_->Save();
	}
	//if (probes_.glvis) { // TODO
	//	std::ofstream omesh("Maxwell1D_RK4.mesh");
	//	omesh.precision(probes_.precision);
	//	mesh_.Print(omesh);
	//	std::ofstream eSol("Maxwell1D_RK4-init-E.gf");
	//	eSol.precision(probes_.precision);
	//	E_.Save(eSol);
	//	std::ofstream hSol("Maxwell1D_RK4-init-H.gf");
	//	hSol.precision(probes_.precision);
	//	H_.Save(hSol);
	//}
}

void Solver::run()
{
	
	setInitialField();

	double time = 0.0;

	maxwellEvol_->SetTime(time);
	odeSolver_->Init(*maxwellEvol_);

	storeInitialVisualizationValues();

	bool done = false;
	int cycle = 0;
	int iter = 0;

	while (!done) {

		odeSolver_->Step(sol_, time, opts_.dt);

		if (probes_.extractDataAtPoints) {
			timeRecord_ = time;
			fieldRecord_ = saveFieldAtPoints(fieldToExtract_);
		}

		done = (time >= opts_.t_final);

		cycle++;

		if (done || cycle % probes_.vis_steps == 0) {
			if (probes_.extractDataAtPoints) {
				timeField_[iter].first = timeRecord_;
				timeField_[iter].second = fieldRecord_;
				iter++;
			}
			if (probes_.paraview) {
				pd_->SetCycle(cycle);
				pd_->SetTime(time);
				pd_->Save();
			}
			//if (probes_.glvis) {
			//	sout_ << "solution\n" << mesh_ << E_ << std::flush; //TODO
			//}
		}
	}
}
}
