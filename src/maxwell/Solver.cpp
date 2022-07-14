#include <fstream>
#include <iostream>
#include <algorithm>

#include "Solver.h"

using namespace mfem;

namespace maxwell {

Solver::Solver(const Model& model, Probes& probes,
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

maxwellEvol_ = std::make_unique<FiniteElementEvolution>(
	fes_.get(), 
	opts_.evolutionOperatorOptions, model_, sources_);

sol_ = Vector(FiniteElementEvolution::numberOfFieldComponents *
	FiniteElementEvolution::numberOfMaxDimensions *
	fes_->GetNDofs());
sol_ = 0.0;



for (int d = X; d <= Z; d++) {
	E_[d].SetSpace(fes_.get());
	E_[d].SetData(sol_.GetData() + d*fes_->GetNDofs());
	H_[d].SetSpace(fes_.get());
	H_[d].SetData(sol_.GetData() + (d+3)*fes_->GetNDofs());
}

initializeSources();

for (int i = 0; i < probes_.getExporterProbes().size(); i++) {
	if (probes_.getExporterProbes().at(i).type == ExporterProbe::Type::Paraview) {
		initializeParaviewData();
		break;
	}
}
for (int i = 0; i < probes_.getExporterProbes().size(); i++) {
	if (probes_.getExporterProbes().at(i).type == ExporterProbe::Type::Glvis) {
		//initializeGLVISData();
		break;
	}
}

if (probes_.getPointsProbes().size()) {
	for (int i = 0; i < probes_.getPointsProbes().size(); i++) {
		elemIds_.resize(probes_.getPointsProbes().size());
		integPointSet_.resize(probes_.getPointsProbes().size());
		auto elemAndIntPointPair = buildElemAndIntegrationPointArrays(probes_.getPointsProbes().at(i).getIntegPointMat());
		elemIds_.at(i) = elemAndIntPointPair.first;
		integPointSet_.at(i) = buildIntegrationPointsSet(elemAndIntPointPair.second);
	}
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

void Solver::initializeSources()
{
	for (int i = 0; i < sources_.getSourcesVector().size(); i++) {
		auto source = sources_.getSourcesVector().at(i);

		std::function<double(const Position&)> f = 0;
		
		switch (model_.getConstMesh().Dimension()) {
		case 1:
			f = std::bind(&Source::evalGaussianFunction1D, &source, std::placeholders::_1);
			break;
		case 2:
			f = std::bind(&Source::evalGaussianFunction2D, &source, std::placeholders::_1);
			break;
		case 3:
			f = std::bind(&Source::evalGaussianFunction3D, &source, std::placeholders::_1);
			break;
		}

		switch (source.getFieldType()) {
		case FieldType::E:
			E_[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
			break;
		case FieldType::H:
			H_[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
			break;
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

const std::pair<Array<int>, Array<IntegrationPoint>> Solver::buildElemAndIntegrationPointArrays(DenseMatrix& physPoints) const
{
	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	std::pair<Array<int>, Array<IntegrationPoint>> res;
	fes_->GetMesh()->FindPoints(physPoints, elemIdArray, integPointArray);
	res = std::make_pair(elemIdArray, integPointArray);
	return res;
}

const std::vector<std::vector<IntegrationPoint>> Solver::buildIntegrationPointsSet(const Array<IntegrationPoint>& ipArray) const
{
	std::vector<IntegrationPoint> aux;
	aux.resize(model_.getConstMesh().Dimension());
	IntegrationPointsSet res;
	res.resize(ipArray.Size(), aux);
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

const std::vector<std::vector<std::array<double, 3>>> Solver::saveFieldAtPointsForAllProbes()
{
	auto maxDir = model_.getConstMesh().Dimension();
	std::vector<FieldFrame> res;
	for (int i = 0; i < probes_.getPointsProbes().size(); i++) {
		FieldFrame aux;
		aux.resize(elemIds_.at(i).Size());
		for (int j = 0; j < elemIds_.at(i).Size(); j++) {
			for (int dir = Direction::X; dir != maxDir; dir++) {
				Direction d = static_cast<Direction>(dir);
				switch (probes_.getPointsProbes().at(i).getFieldType()) {
				case FieldType::E:
					aux[j][probes_.getPointsProbes().at(i).getDirection()] =
						E_[probes_.getPointsProbes().at(i).getDirection()].GetValue(elemIds_.at(i)[j], integPointSet_.at(i).at(j)[d]);
					break;
				case FieldType::H:
					aux[j][probes_.getPointsProbes().at(i).getDirection()] =
						H_[probes_.getPointsProbes().at(i).getDirection()].GetValue(elemIds_.at(i)[j], integPointSet_.at(i).at(j)[d]);
					break;
				}
			}
		}
		res.push_back(aux);
	}
	return res;
}

void Solver::initializeParaviewData()
{
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

void Solver::storeInitialVisualizationValues()
{
	for (int i = 0; i < probes_.getExporterProbes().size(); i++) {
		if (probes_.getExporterProbes().at(i).type == ExporterProbe::Type::Paraview) {
			pd_->SetCycle(0);
			pd_->SetTime(0.0);
			pd_->Save();
			break;
		}
	}
}

void Solver::run()
{

	double time = 0.0;

	maxwellEvol_->SetTime(time);
	odeSolver_->Init(*maxwellEvol_);

	storeInitialVisualizationValues();

	bool done = false;
	int cycle = 0;
	
	if (probes_.getPointsProbes().size()) {
		timeRecord_ = time;
		fieldRecord_ = saveFieldAtPointsForAllProbes();
		for (int i = 0; i < probes_.getPointsProbes().size(); i++) {
			probes_.getPointsProbes().at(i).getFieldMovie().emplace(timeRecord_, fieldRecord_.at(i));
		}
	}

	while (!done) {

		odeSolver_->Step(sol_, time, opts_.dt);

		if (abs(time - opts_.t_final) < 1e-6) {
			done = true;
		}

		cycle++;

		if (done || cycle % probes_.vis_steps == 0) {
			if (probes_.getPointsProbes().size()) {
				fieldRecord_ = saveFieldAtPointsForAllProbes();
				for (int i = 0; i < probes_.getPointsProbes().size(); i++) {
					timeRecord_ = time;
					probes_.getPointsProbes().at(i).getFieldMovie().emplace(timeRecord_, fieldRecord_.at(i));
				}
			}
			for (int i = 0; i < probes_.getExporterProbes().size(); i++) {
				if (probes_.getExporterProbes().at(i).type == ExporterProbe::Type::Paraview) {
				pd_->SetCycle(cycle);
				pd_->SetTime(time);
				pd_->Save();
				break;
				}
			}
			//if (probes_.glvis) {
			//	sout_ << "solution\n" << mesh_ << E_ << std::flush; //TODO
			//}
		}
	}
}
}
