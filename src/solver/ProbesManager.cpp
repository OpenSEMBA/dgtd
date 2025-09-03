#include "ProbesManager.h"
#include "math/PhysicalConstants.h"

namespace maxwell {

using namespace mfem;

std::string getFieldPolString(const FieldType& ft, const Direction& d)
{
	switch(ft){
		case E:
		switch(d){
			case X:
				return "Ex";
			case Y:
				return "Ey";
			case Z:
				return "Ez";
			default:
				throw std::runtime_error("Incorrect direction in getFieldPolString.");
		}
		case H:
		switch(d){
			case X:
				return "Hx";
			case Y:
				return "Hy";
			case Z:
				return "Hz";
			default:
				throw std::runtime_error("Incorrect direction in getFieldPolString.");
		}
		default:
			throw std::runtime_error("Incorrect fieldtype in getFieldPolString.");
	}
}

ParaViewDataCollection ProbesManager::buildParaviewDataCollectionInfo(const ExporterProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const
{
	fes_.ExchangeFaceNbrData();
	fes_.GetParMesh()->ExchangeFaceNbrData();
	ParaViewDataCollection pd{ p.name, fes_.GetParMesh()};
	pd.SetPrefixPath("ParaView");

	pd.RegisterField("E", &fields.get(E));
	pd.RegisterField("H", &fields.get(H));
	
	bool highOrder = false;
	auto geomElemOrder = fes_.GetMesh()->GetElementTransformation(0)->Order();
	auto fecorder = fes_.FEColl()->GetOrder();
	geomElemOrder > 1 || fecorder > 1 ? highOrder = true : highOrder = false;
	pd.SetHighOrderOutput(highOrder);
	pd.SetLevelsOfDetail(fecorder);
	
	pd.SetDataFormat(VTKFormat::BINARY);

	return pd;
}

ProbesManager::ProbesManager(Probes pIn, mfem::ParFiniteElementSpace& fes, Fields<ParFiniteElementSpace, ParGridFunction>& fields, const SolverOptions& opts) :
	probes{ pIn },
	fes_{ fes }
{
	for (const auto& p: probes.exporterProbes) {
		exporterProbesCollection_.emplace(&p, buildParaviewDataCollectionInfo(p, fields));
	}

	for (const auto& p : probes.pointProbes) {
		pointProbesCollection_.emplace(&p, buildPointProbeCollectionInfo(p, fields));
	}

	for (const auto& p : probes.fieldProbes) {
		fieldProbesCollection_.emplace(&p, buildFieldProbeCollectionInfo(p, fields));
	}

	for (const auto& p : probes.nearFieldProbes) {
		auto dgfec{ dynamic_cast<const DG_FECollection*>(fes_.FEColl()) };
		if (!dgfec)
		{
			throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
		}
		nearFieldReqs_.emplace(&p, std::make_unique<NearFieldReqs>(NearFieldReqs(p, dgfec, fes_, fields)));
		nearFieldProbesCollection_.emplace(&p, buildNearFieldDataCollectionInfo(p, fields));
	}

	for (const auto& p: probes.domainSnapshotProbes) {
		domainSnapshotProbesCollection_.emplace(&p, buildDomainSnapshotDataCollection(p, fields));
	}
	
	finalTime_ = opts.finalTime;
	fields_ = &fields;
}

const FieldProbe& ProbesManager::getFieldProbe(const std::size_t i) const
{
	assert(i < probes.fieldProbes.size());
	return probes.fieldProbes[i];
}

const PointProbe& ProbesManager::getPointProbe(const std::size_t i) const
{
	assert(i < probes.pointProbes.size());
	return probes.pointProbes[i];
}

const ParGridFunction& getFieldView(const FieldProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields)
{
	switch (p.getFieldType()) {
	case FieldType::E:
		return fields.get(E, p.getDirection());
	case FieldType::H:
		return fields.get(H, p.getDirection());
	default:
		throw std::runtime_error("Invalid field type.");
	}
}

DenseMatrix pointVectorToDenseMatrixColumnVector(const Point& p)
{
	DenseMatrix r{(int)p.size(), 1 };
	for (auto i{ 0 }; i < p.size(); ++i) {
		r(i, 0) = p[i];
	}
	return r;
}

ProbesManager::PointProbeCollection
ProbesManager::buildPointProbeCollectionInfo(const PointProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const
{
	
	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	auto pointMatrix{ pointVectorToDenseMatrixColumnVector(p.getPoint()) };
	fes_.GetParMesh()->FindPoints(pointMatrix, elemIdArray, integPointArray);
	assert(elemIdArray.Size() == 1);
	assert(integPointArray.Size() == 1);
	FESPoint fesPoints { elemIdArray[0], integPointArray[0] };

	return { 
		fesPoints, 
		fields.get(E, X),
		fields.get(E, Y),
		fields.get(E, Z),
		fields.get(H, X),
		fields.get(H, Y),
		fields.get(H, Z)
	};
}

void ProbesManager::initPointFieldProbeExport()
{
	for (const auto& p : probes.pointProbes) {
		if(p.write){
			if (Mpi::WorldRank() == 0){
				std::ofstream myfile;
				std::string path("SimulationExports/" + caseName_ + "/" + "PointProbe" + std::to_string(p.getProbeID()) + ".dat");
				std::vector<double> position = std::vector<double>({0.0, 0.0, 0.0});
				for (auto i = 0; i < p.getPoint().size(); i++){
					position[i] = p.getPoint()[i];
				}
				myfile.open(path, std::ios::app);
				if (myfile.is_open()) {
					myfile << "PointProbe ID " << std::to_string(p.getProbeID()) << "\n";
					myfile << "Spatial Position (X, Y, Z) \n";
					myfile << std::scientific << std::setprecision(5);
					myfile << std::to_string(position[0]) + " " + std::to_string(position[1]) + " " + std::to_string(position[2]) << "\n";
					myfile << "Time (s) // Ex // Ey // Ez // Hx // Hy // Hz \n";
				}
				myfile.close();
			}
		}
	}

	for (const auto& p : probes.fieldProbes) {
		if(p.write){
			if (Mpi::WorldRank() == 0){
				std::ofstream myfile;
				std::string path("SimulationExports/" + caseName_ + "/" + "FieldProbe" + std::to_string(p.getProbeID()) + ".dat");
				std::vector<double> position = std::vector<double>({0.0, 0.0, 0.0});
				auto fieldpol = getFieldPolString(p.getFieldType(), p.getDirection());
				for (auto i = 0; i < p.getPoint().size(); i++){
					position[i] = p.getPoint()[i];
				}
				myfile.open(path, std::ios::app);
				if (myfile.is_open()) {
					myfile << "FieldProbe ID " << std::to_string(p.getProbeID()) << "\n";
					myfile << "Spatial Position (X, Y, Z) \n";
					myfile << std::scientific << std::setprecision(5);
					myfile << std::to_string(position[0]) + " " + std::to_string(position[1]) + " " + std::to_string(position[2]) << "\n";
					myfile << "Time (s) // " + fieldpol + "\n";
				}
				myfile.close();
			}
		}
	}
}

ProbesManager::FieldProbeCollection
ProbesManager::buildFieldProbeCollectionInfo(const FieldProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const
{

	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	auto pointMatrix{ pointVectorToDenseMatrixColumnVector(p.getPoint()) };
	fes_.GetParMesh()->FindPoints(pointMatrix, elemIdArray, integPointArray);
	assert(elemIdArray.Size() == 1);
	assert(integPointArray.Size() == 1);
	FESPoint fesPoints{ elemIdArray[0], integPointArray[0] };

	return {
		fesPoints,
		getFieldView(p, fields)
	};
}

void isDGCollection(const FiniteElementSpace& fes)
{
	if (!dynamic_cast<const DG_FECollection*>(fes.FEColl()))
	{
		throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
	}
}

DataCollection ProbesManager::buildNearFieldDataCollectionInfo(
	const NearFieldProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& gFields) const
{

	isDGCollection(fes_);

	DataCollection res{ p.name, nearFieldReqs_.at(&p)->getSubMesh() };
	res.SetPrefixPath("NearToFarFieldExports/" + p.name + "/rank" + std::to_string(Mpi::WorldRank()));
	res.RegisterField("Ex.gf", &nearFieldReqs_.at(&p)->getField(E, X));
	res.RegisterField("Ey.gf", &nearFieldReqs_.at(&p)->getField(E, Y));
	res.RegisterField("Ez.gf", &nearFieldReqs_.at(&p)->getField(E, Z));
	res.RegisterField("Hx.gf", &nearFieldReqs_.at(&p)->getField(H, X));
	res.RegisterField("Hy.gf", &nearFieldReqs_.at(&p)->getField(H, Y));
	res.RegisterField("Hz.gf", &nearFieldReqs_.at(&p)->getField(H, Z));

	return res;

}

DataCollection ProbesManager::buildDomainSnapshotDataCollection(const DomainSnapshotProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields)
{

	isDGCollection(fes_);

	DataCollection res{ p.name, fes_.GetParMesh() };
	res.SetFormat(res.PARALLEL_FORMAT);
	res.SetPrefixPath("DomainSnapshotExports/" + p.name + "/rank" + std::to_string(Mpi::WorldRank()));
	res.RegisterField("Ex.gf", &fields.get(E, X));
	res.RegisterField("Ey.gf", &fields.get(E, Y));
	res.RegisterField("Ez.gf", &fields.get(E, Z));
	res.RegisterField("Hx.gf", &fields.get(H, X));
	res.RegisterField("Hy.gf", &fields.get(H, Y));
	res.RegisterField("Hz.gf", &fields.get(H, Z));

	const auto& ex = &fields.get(E, X);

	return res;
}


void ProbesManager::updateProbe(ExporterProbe& p, Time time)
{
	if (std::abs(time - finalTime_) >= 1e-3){
		if (cycle_ % p.visSteps != 0) {
			return;
		}
	}

	auto it{ exporterProbesCollection_.find(&p) };
	assert(it != exporterProbesCollection_.end());
	auto& pd{ it->second };

	pd.SetCycle(cycle_);
	pd.SetTime(time);

	pd.Save();
}

void ProbesManager::updateProbe(FieldProbe& p, Time time)
{
	const auto& it{ fieldProbesCollection_.find(&p) };
	assert(it != fieldProbesCollection_.end());
	const auto& pC{ it->second };
	if (pC.fesPoint.elementId != -2){

		real_t gf_value = pC.field.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);

		p.addFieldToMovies(
			time, 
			gf_value
		);

		if(p.write){
			std::ofstream myfile;
			std::string path("SimulationExports/" + caseName_ + "/" + "FieldProbe" + std::to_string(p.getProbeID()) + ".dat");
			myfile.open(path, std::ios::app);
			if (myfile.is_open()) {
				myfile << std::scientific << std::setprecision(5);
				myfile << time / physicalConstants::speedOfLight_SI << " " << gf_value << "\n";
			}
		myfile.close();
		}
	}
}

void ProbesManager::updateProbe(PointProbe& p, Time time)
{
	const auto& it{ pointProbesCollection_.find(&p) };
	assert(it != pointProbesCollection_.end());
	const auto& pC{ it->second };
	if (pC.fesPoint.elementId != -2){
		FieldsForMovie f4FP;
		{
			f4FP.Ex = pC.field_Ex.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
			f4FP.Ey = pC.field_Ey.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
			f4FP.Ez = pC.field_Ez.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
			f4FP.Hx = pC.field_Hx.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
			f4FP.Hy = pC.field_Hy.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
			f4FP.Hz = pC.field_Hz.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
		}
		p.addFieldsToMovies(
			time,
			f4FP
		);
		if(p.write){
			std::ofstream myfile;
			std::string path("SimulationExports/" + caseName_ + "/" + "PointProbe" + std::to_string(p.getProbeID()) + ".dat");
			myfile.open(path, std::ios::app);
			if (myfile.is_open()) {
				myfile << std::scientific << std::setprecision(5);
				myfile << time / physicalConstants::speedOfLight_SI << 
				" " << f4FP.Ex << " " << f4FP.Ey << " " << f4FP.Ez <<
				" " << f4FP.Hx << " " << f4FP.Hy << " " << f4FP.Hz << "\n";
			}
			myfile.close();
		}
	}
}

Fields<ParFiniteElementSpace, ParGridFunction> buildFieldsForProbe(const Fields<ParFiniteElementSpace, ParGridFunction>& src, ParFiniteElementSpace& fes)
{
	Fields<ParFiniteElementSpace, ParGridFunction> res(fes);
	for (auto f : { E, H }) {
		for (auto& d : { X, Y, Z }) {
			TransferMap tm(src.get(f, d), res.get(f, d));
			tm.   Transfer(src.get(f, d), res.get(f, d));
		}
	}
	return res;
}

void ProbesManager::updateProbe(NearFieldProbe& p, Time time)
{
	if (std::abs(time - finalTime_) >= 1e-3) {
		if (cycle_ % p.expSteps != 0) {
			return;
		}
	}

	auto it{ nearFieldProbesCollection_.find(&p) };
	assert(it != nearFieldProbesCollection_.end());
	auto& dc{ it->second };

	nearFieldReqs_.at(&p)->updateFields();

	dc.SetCycle(cycle_);
	dc.SetTime(time);
	dc.Save(); //Ideally we wouldn't save the mesh every step but MFEM DataCollection does it. We should redo the method or delete the extra meshes later (suboptimal).

	std::string mesh_path{ dc.GetPrefixPath() + dc.GetCollectionName() + "/mesh" }; //We do want to save the mesh at the base directory, though ideally we'd only do this once. (WIP)
	auto mesh{ dc.GetMesh() };
	auto elemOrder = fes_.GetMesh()->GetElementTransformation(0)->Order();
	mesh->SetCurvature(elemOrder);
	mesh->Save(dc.GetPrefixPath() + "/mesh");

	std::string dir_name = dc.GetPrefixPath() + dc.GetCollectionName() + "_" + to_padded_string(dc.GetCycle(), 6) + "/time.txt";
	std::ofstream file;
	file.open(dir_name);
	file << time;
	file.close();
}

void ProbesManager::updateProbe(DomainSnapshotProbe& p, Time time)
{
	if (std::abs(time - finalTime_) >= 1e-3) {
		if (cycle_ % p.expSteps != 0) {
			return;
		}
	}

	auto it{ domainSnapshotProbesCollection_.find(&p) };
	assert(it != domainSnapshotProbesCollection_.end());
	auto& dc{ it->second };

	dc.SetCycle(cycle_);
	dc.SetTime(time);
	dc.SetFormat(1); // 0 - Serial Format / 1 - Parallel Format
	dc.Save();

	std::string dir_name = dc.GetPrefixPath() + dc.GetCollectionName() + "_" + to_padded_string(dc.GetCycle(), 6) + "/time.txt";
	std::ofstream file;
	file.open(dir_name);
	file << time;
	file.close();
	
}

void ProbesManager::updateProbes(Time t)
{
	for (auto& p : probes.exporterProbes) {
		updateProbe(p, t);
	}
	
	for (auto& p : probes.fieldProbes) {
		updateProbe(p, t);
	}

	for (auto& p : probes.pointProbes) {
		updateProbe(p, t);
	}
	
	for (auto& p : probes.nearFieldProbes) {
		updateProbe(p, t);
	}

	for (auto& p : probes.domainSnapshotProbes){
		updateProbe(p, t);
	}

	cycle_++;
}

Array<int> buildSurfaceMarker(const std::vector<int>& tags, const ParFiniteElementSpace& fes)
{
	Array<int> res(fes.GetMesh()->bdr_attributes.Max());
	res = 0;
	for (const auto& t : tags) {
		res[t - 1] = 1;
	}
	return res;
}

void NearFieldReqs::updateFields()
{
	tMaps_.transferFields(gFields_, fields_);
}

NearFieldReqs::NearFieldReqs(
	const NearFieldProbe& p, const DG_FECollection* fec, ParFiniteElementSpace& fes, Fields<ParFiniteElementSpace, ParGridFunction>& global) :
	ntff_smsh_{ NearToFarFieldSubMesher(*fes.GetMesh(), fes, buildSurfaceMarker(p.tags, fes)) },
	sfes_{ std::make_unique<FiniteElementSpace>(ntff_smsh_.getSubMesh(), fec) },
	fields_{ Fields<FiniteElementSpace, GridFunction>(*sfes_) },
	gFields_{ global },
	tMaps_{ TransferMaps(gFields_, fields_) }
{
	updateFields();
}

}
