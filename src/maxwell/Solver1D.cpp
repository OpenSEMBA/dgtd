#include <fstream>
#include <iostream>
#include <algorithm>

#include "Solver1D.h"

using namespace mfem;

namespace maxwell {

	Solver1D::Solver1D(const Model& model, const Probes& probes,
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

	maxwellEvol_ = std::make_unique<FiniteElementEvolutionNoCond>(fes_.get(), opts_.evolutionOperatorOptions);

	sol_ = Vector(FiniteElementEvolutionNoCond::numberOfFieldComponents * mesh_.Dimension() * fes_->GetNDofs());
	sol_ = 0.0;

	for (int d = X; d <= Z; d++) {
		E_[d].SetSpace(fes_.get());
		E_[d].SetData(sol_.GetData() + d*fes_->GetNDofs());
		H_[d].SetSpace(fes_.get());
		H_[d].SetData(sol_.GetData() + (d+3)*fes_->GetNDofs());
	}


	if (probes_.paraview) {
		initializeParaviewData();
	}
	if (probes_.glvis) {
		initializeGLVISData();
	}
	if (probes_.extractDataAtPoint) {
		integPoint_ = setIntegrationPoint(probes_.integPoint);
		fieldToExtract_ = probes_.fieldToExtract;
	}
}

void Solver1D::checkOptionsAreValid(const Options& opts, const Mesh& mesh)
{
	if (mesh.Dimension() != 1) {
		throw std::exception("Incorrect Dimension for mesh");
	}
	if ((opts.order < 0) ||
		(opts.t_final < 0) ||
		(opts.dt < 0)) {
		throw std::exception("Incorrect parameters in Options");
	}
}


void Solver1D::setInitialField(const FieldType& ft, std::function<double(const Position&)> f, const Direction& d)
{
	switch (ft) {
	case FieldType::E:
		E_[d].ProjectCoefficient(FunctionCoefficient(f));
		return;
	case FieldType::H:
		H_[d].ProjectCoefficient(FunctionCoefficient(f));
		return;
	}
}

const GridFunction& Solver1D::getFieldInDirection(const FieldType& ft, const Direction& d) const
{
	switch (ft) {
	case FieldType::E:
		return E_[d];
	case FieldType::H:
		return H_[d];
	}
}

const Vector& Solver1D::getMaterialProperties(const Material& mat) const
{
	return Vector({mat.getPermittivity(), mat.getPermeability(), mat.getImpedance(), mat.getConductance()});
}

const int Solver1D::getElementIndexForPosition(const IntegrationPoint& ip) const
{
	Vector meshBoundingMin, meshBoundingMax;
	Mesh meshc = Mesh(mesh_, true);
	meshc.GetBoundingBox(meshBoundingMin, meshBoundingMax);
	int res = (int)std::floor(ip.x / ((meshBoundingMax[0] - meshBoundingMin[0]) / meshc.GetNE())) - 1;
	return res;
}

const Array<double> Solver1D::getVertexPositionInPhysicalCoords(const Array<int>& elementVertex) const
{
	Vector meshBoundingMin, meshBoundingMax;
	Mesh meshc = Mesh(mesh_, true);
	meshc.GetBoundingBox(meshBoundingMin, meshBoundingMax);
	Array<double> res(elementVertex.Size());
	const double vertexTop = meshc.GetNV() - 1;
	for (int i = 0; i < elementVertex.Size(); i++) {
		res[i] = (1.0 - ((vertexTop - elementVertex[i]) / vertexTop)) / (meshBoundingMax[0] - meshBoundingMin[0]);
	}
	return res;


}

const IntegrationPoint Solver1D::getRelativePositionInElement(const int& elementN, const IntegrationPoint& ip) const
{
	Array<int> aux;
	mesh_.GetElementVertices(elementN, aux);
	Array<double> auxPos = getVertexPositionInPhysicalCoords(aux);

	IntegrationPoint res;
	res.Set1w(((ip.x - auxPos[0]) / (auxPos[1] - auxPos[0])), 0.0);
	return res;
}

//const double Solver1D::saveFieldAtPoint(const IntegrationPoint& ip, const FieldType& ft) const //TODO
//{
//	/*switch (ft) {*/ //TODO EXTEND 3D
//	//case FieldType::E:
//	//	return E_.GetValue(getElementIndexForPosition(ip),
//	//		getRelativePositionInElement(getElementIndexForPosition(ip), ip));
//	//case FieldType::H:
//	//	return H_.GetValue(getElementIndexForPosition(ip),
//	//		getRelativePositionInElement(getElementIndexForPosition(ip), ip));
//	//}
//}

const IntegrationPoint Solver1D::setIntegrationPoint(const IntegrationPoint& ip) const
{
	IntegrationPoint res;
	res.Set1w(ip.x, ip.weight);
	return res;
}

void Solver1D::initializeParaviewData()
{
	pd_ = NULL;
	pd_ = std::make_unique<ParaViewDataCollection>("MaxwellView1D", &mesh_);
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

void Solver1D::initializeGLVISData()
{
	char vishost[] = "localhost";
	int  visport = 19916;
	sout_.open(vishost, visport);
	sout_.precision(probes_.precision);
	//sout_ << "solution\n" << mesh_ << E_; //TODO
	sout_ << "pause\n";
	sout_ << std::flush;
	std::cout << "GLVis visualization paused."
		<< " Press space (in the GLVis window) to resume it.\n";
}

void Solver1D::storeInitialVisualizationValues()
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

void Solver1D::run()
{
	Vector vals(2 * fes_.get()->GetVSize());

	double time = 0.0;

	maxwellEvol_->SetTime(time);
	odeSolver_->Init(*maxwellEvol_);

	storeInitialVisualizationValues();

	if (probes_.extractDataAtPoint) {
		timeRecord_.SetSize(std::ceil(opts_.t_final / opts_.dt));
		fieldRecord_.SetSize(std::ceil(opts_.t_final / opts_.dt));
	}

	bool done = false;
	int cycle = 0;

	while (!done) {
		odeSolver_->Step(sol_, time, opts_.dt);

		if (probes_.extractDataAtPoint) {
			timeRecord_[cycle] = time;
			//fieldRecord_[cycle] = saveFieldAtPoint(integPoint_, fieldToExtract_);
		}

		done = (time >= opts_.t_final);

		cycle++;

		if (done || cycle % probes_.vis_steps == 0) {
			if (probes_.extractDataAtPoint) {

				timeField_.SetSize(std::ceil(opts_.t_final / opts_.dt) * 2);

				for (int i = 0; i < timeField_.Size()/2; i++) {

					timeField_.Elem(i) = timeRecord_.Elem(i);
					//timeField_.Elem(i + timeField_.Size() / 2) = fieldRecord_.Elem(i);

				}
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
