#include "Solver1D.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace maxwell {

Solver1D::Solver1D(const Options& opts, const Mesh& mesh)
{
	checkOptionsAreValid(opts, mesh);

	mesh_ = mfem::Mesh(mesh, true);
	opts_ = opts;

	fec_ = std::make_unique<DG_FECollection>(
		opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
	fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	odeSolver_ = std::make_unique<RK4Solver>();

	maxwellEvol_ = std::make_unique<FE_Evolution>(fes_.get(), opts.evolutionOperatorOptions);

	sol_ = Vector(FE_Evolution::numberOfFieldComponents * fes_->GetNDofs());
	sol_ = 0.0;

	E_.SetSpace(fes_.get());
	E_.SetData(sol_.GetData());
	H_.SetSpace(fes_.get());
	H_.SetData(sol_.GetData() + fes_->GetNDofs());

	if (opts_.paraview) {
		initializeParaviewData();
	}
	if (opts_.glvis) {
		initializeGLVISData();
	}
	if (opts_.extractDataAtPoint) {
		integPoint_ = setIntegrationPoint(opts_.integPoint);
		fieldToExtract_ = opts_.fieldToExtract;
	}
}

void Solver1D::checkOptionsAreValid(const Options& opts, const Mesh& mesh)
{
	if (mesh.Dimension() != 1) {
		throw std::exception("Incorrect Dimension for mesh");
	}
	if ((opts.order < 0) ||
		(opts.t_final < 0) ||
		(opts.dt < 0) ||
		(opts.vis_steps < 1) ||
		(opts.precision < 1)) {
		throw std::exception("Incorrect parameters in Options");
	}
}

mfem::Array<int> Solver1D::buildEssentialTrueDOF()
{
	Array<int> ess_tdof_list;
	if (mesh_.bdr_attributes.Size())
	{
		Array<int> ess_bdr(mesh_.bdr_attributes.Max());
		ess_bdr = 1;
		fes_.get()->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
	}
	return ess_tdof_list;
}

void Solver1D::setInitialField(const FieldType& ft, std::function<double(const Position&)> f)
{
	switch (ft) {
	case FieldType::Electric:
		E_.ProjectCoefficient(FunctionCoefficient(f));
		return;
	case FieldType::Magnetic:
		H_.ProjectCoefficient(FunctionCoefficient(f));
		return;
	}
}

const GridFunction& Solver1D::getField(const FieldType& ft) const
{
	switch (ft) {
	case FieldType::Electric:
		return E_;
	case FieldType::Magnetic:
		return H_;
	}
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

const double Solver1D::saveFieldAtPoint(const IntegrationPoint& ip, const FieldType& ft) const
{
	switch (ft) {
	case FieldType::Electric:
		return E_.GetValue(getElementIndexForPosition(ip),
			getRelativePositionInElement(getElementIndexForPosition(ip), ip));
	case FieldType::Magnetic:
		return H_.GetValue(getElementIndexForPosition(ip),
			getRelativePositionInElement(getElementIndexForPosition(ip), ip));
	}
}

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
	pd_->RegisterField("E", &E_);
	pd_->RegisterField("H", &H_);
	pd_->SetLevelsOfDetail(opts_.order);
	pd_->SetDataFormat(VTKFormat::BINARY);
	opts_.order > 0 ? pd_->SetHighOrderOutput(true) : pd_->SetHighOrderOutput(false);
}

void Solver1D::initializeGLVISData()
{
	char vishost[] = "localhost";
	int  visport = 19916;
	sout_.open(vishost, visport);
	sout_.precision(opts_.precision);
	sout_ << "solution\n" << mesh_ << E_;
	sout_ << "pause\n";
	sout_ << std::flush;
	std::cout << "GLVis visualization paused."
		<< " Press space (in the GLVis window) to resume it.\n";
}

void Solver1D::storeInitialVisualizationValues()
{
	if (opts_.paraview) {
		pd_->SetCycle(0);
		pd_->SetTime(0.0);
		pd_->Save();
	}
	if (opts_.glvis) {
		std::ofstream omesh("Maxwell1D_RK4.mesh");
		omesh.precision(opts_.precision);
		mesh_.Print(omesh);
		std::ofstream eSol("Maxwell1D_RK4-init-E.gf");
		eSol.precision(opts_.precision);
		E_.Save(eSol);
		std::ofstream hSol("Maxwell1D_RK4-init-H.gf");
		hSol.precision(opts_.precision);
		H_.Save(hSol);
	}
}

void Solver1D::run()
{
	Vector vals(2 * fes_.get()->GetVSize());

	double time = 0.0;

	maxwellEvol_->SetTime(time);
	odeSolver_->Init(*maxwellEvol_);

	storeInitialVisualizationValues();

	if (opts_.extractDataAtPoint) {
		timeRecord_.SetSize(std::ceil(opts_.t_final / opts_.dt));
		fieldRecord_.SetSize(std::ceil(opts_.t_final / opts_.dt));
	}

	bool done = false;
	int cycle = 0;

	while (!done) {
		odeSolver_->Step(sol_, time, opts_.dt);

		if (opts_.extractDataAtPoint) {
			timeRecord_[cycle] = time;
			fieldRecord_[cycle] = saveFieldAtPoint(integPoint_, fieldToExtract_);
		}

		done = (time >= opts_.t_final);

		cycle++;

		if (done || cycle % opts_.vis_steps == 0) {
			if (opts_.extractDataAtPoint) {

				timeField_.SetSize(std::ceil(opts_.t_final / opts_.dt) * 2);

				for (int i = 0; i < timeField_.Size()/2; i++) {

					timeField_.Elem(i) = timeRecord_.Elem(i);
					timeField_.Elem(i + timeField_.Size() / 2) = fieldRecord_.Elem(i);

				}
			}
			if (opts_.paraview) {
				pd_->SetCycle(cycle);
				pd_->SetTime(time);
				pd_->Save();
			}
			if (opts_.glvis) {
				sout_ << "solution\n" << mesh_ << E_ << std::flush;
			}
		}
	}
}
}
