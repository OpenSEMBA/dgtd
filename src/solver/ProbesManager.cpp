#include "ProbesManager.h"

namespace maxwell {

using namespace mfem;

ParaViewDataCollection ProbesManager::buildParaviewDataCollectionInfo(const ExporterProbe& p, Fields& fields) const
{
	ParaViewDataCollection pd{ p.name, fes_.GetMesh()};
	pd.SetPrefixPath("ParaView");
	
	pd.RegisterField("Ex", &fields.get(E,X));
	pd.RegisterField("Ey", &fields.get(E,Y));
	pd.RegisterField("Ez", &fields.get(E,Z));
	pd.RegisterField("Hx", &fields.get(H,X));
	pd.RegisterField("Hy", &fields.get(H,Y));
	pd.RegisterField("Hz", &fields.get(H,Z));
	
	const auto order{ fes_.GetMaxElementOrder() };
	pd.SetLevelsOfDetail(3);
	pd.SetHighOrderOutput(true);
	
	pd.SetDataFormat(VTKFormat::BINARY);

	return pd;
}

ProbesManager::ProbesManager(Probes pIn, mfem::FiniteElementSpace& fes, Fields& fields, const SolverOptions& opts) :
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

	for (const auto& p : probes.nearToFarFieldProbes) {
		auto dgfec{ dynamic_cast<const DG_FECollection*>(fes_.FEColl()) };
		if (!dgfec)
		{
			throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
		}
		nearToFarFieldReqs_.emplace(&p, std::make_unique<NearToFarFieldReqs>(NearToFarFieldReqs(p, dgfec, fes_, fields)));
		nearToFarFieldProbesCollection_.emplace(&p, buildNearToFarFieldDataCollectionInfo(p, fields));
	}

	finalTime_ = opts.finalTime;
}

const PointProbe& ProbesManager::getPointProbe(const std::size_t i) const
{
	assert(i < probes.pointProbes.size());
	return probes.pointProbes[i];
}

const FieldProbe& ProbesManager::getFieldProbe(const std::size_t i) const
{
	assert(i < probes.fieldProbes.size());
	return probes.fieldProbes[i];
}

const GridFunction& getFieldView(const PointProbe& p, Fields& fields)
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
ProbesManager::buildPointProbeCollectionInfo(const PointProbe& p, Fields& fields) const
{
	
	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	auto pointMatrix{ pointVectorToDenseMatrixColumnVector(p.getPoint()) };
	fes_.GetMesh()->FindPoints(pointMatrix, elemIdArray, integPointArray);
	assert(elemIdArray.Size() == 1);
	assert(integPointArray.Size() == 1);
	FESPoint fesPoints { elemIdArray[0], integPointArray[0] };
	
	return { 
		fesPoints, 
		getFieldView(p, fields)
	};
}

ProbesManager::FieldProbeCollection
ProbesManager::buildFieldProbeCollectionInfo(const FieldProbe& p, Fields& fields) const
{

	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	auto pointMatrix{ pointVectorToDenseMatrixColumnVector(p.getPoint()) };
	fes_.GetMesh()->FindPoints(pointMatrix, elemIdArray, integPointArray);
	assert(elemIdArray.Size() == 1);
	assert(integPointArray.Size() == 1);
	FESPoint fesPoints{ elemIdArray[0], integPointArray[0] };

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

DataCollection ProbesManager::buildNearToFarFieldDataCollectionInfo(
	const NearToFarFieldProbe& p, Fields& gFields) const
{
	if (!dynamic_cast<const DG_FECollection*>(fes_.FEColl()))
	{
		throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
	}

	DataCollection res{ p.name, nearToFarFieldReqs_.at(&p)->getSubMesh()};
	res.SetPrefixPath("NearToFarFieldExports/" + p.name);
	res.RegisterField("Ex.gf", &nearToFarFieldReqs_.at(&p)->getField(E, X));
	res.RegisterField("Ey.gf", &nearToFarFieldReqs_.at(&p)->getField(E, Y));
	res.RegisterField("Ez.gf", &nearToFarFieldReqs_.at(&p)->getField(E, Z));
	res.RegisterField("Hx.gf", &nearToFarFieldReqs_.at(&p)->getField(H, X));
	res.RegisterField("Hy.gf", &nearToFarFieldReqs_.at(&p)->getField(H, Y));
	res.RegisterField("Hz.gf", &nearToFarFieldReqs_.at(&p)->getField(H, Z));

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

void ProbesManager::updateProbe(PointProbe& p, Time time)
{
	const auto& it{ pointProbesCollection_.find(&p) };
	assert(it != pointProbesCollection_.end());
	const auto& pC{ it->second };
	p.addFieldToMovies(
		time, 
		pC.field.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP)
	);
}

void ProbesManager::updateProbe(FieldProbe& p, Time time)
{
	const auto& it{ fieldProbesCollection_.find(&p) };
	assert(it != fieldProbesCollection_.end());
	const auto& pC{ it->second };
	FieldsForFP f4FP;
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
}

Fields buildFieldsForProbe(const Fields& src, FiniteElementSpace& fes)
{
	Fields res(fes);
	for (auto f : { E, H }) {
		for (auto& d : { X, Y, Z }) {
			TransferMap tm(src.get(f, d), res.get(f, d));
			tm.   Transfer(src.get(f, d), res.get(f, d));
		}
	}
	return res;
}

void ProbesManager::updateProbe(NearToFarFieldProbe& p, Time time)
{
	if (abs(time - finalTime_) >= 1e-3) {
		if (cycle_ % p.steps != 0) {
			return;
		}
	}

	auto it{ nearToFarFieldProbesCollection_.find(&p) };
	assert(it != nearToFarFieldProbesCollection_.end());
	auto& dc{ it->second };

	nearToFarFieldReqs_.at(&p)->updateFields();

	dc.SetCycle(cycle_);
	dc.SetTime(time);
	dc.Save(); //Ideally we wouldn't save the mesh every step but MFEM DataCollection does it. We should redo the method or delete the extra meshes later (suboptimal).

	std::string mesh_path{ dc.GetPrefixPath() + dc.GetCollectionName() + "/mesh" }; //We do want to save the mesh at the base directory, though ideally we'd only do this once. (WIP)
	auto mesh{ dc.GetMesh() };
	mesh->Save(dc.GetPrefixPath() + "/mesh");

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
	
	for (auto& p : probes.pointProbes) {
		updateProbe(p, t);
	}

	for (auto& p : probes.fieldProbes) {
		updateProbe(p, t);
	}

	for (auto& p : probes.nearToFarFieldProbes) {
		updateProbe(p, t);
	}

	cycle_++;
}

Array<int> buildSurfaceMarker(const std::vector<int>& tags, const FiniteElementSpace& fes)
{
	Array<int> res(fes.GetMesh()->bdr_attributes.Max());
	res = 0;
	for (const auto& t : tags) {
		res[t - 1] = 1;
	}
	return res;
}

void NearToFarFieldReqs::assignGlobalFieldsReferences(Fields& global)
{
	gFields_.get(E, X) = global.get(E, X);
	gFields_.get(E, Y) = global.get(E, Y);
	gFields_.get(E, Z) = global.get(E, Z);
	gFields_.get(H, X) = global.get(H, X);
	gFields_.get(H, Y) = global.get(H, Y);
	gFields_.get(H, Z) = global.get(H, Z);
}

void NearToFarFieldReqs::updateFields()
{
	tMaps_.transferFields(gFields_, fields_);
}

void TransferMaps::transferFields(const Fields& src, Fields& dst)
{
	tMapEx.Transfer(src.get(E, X), dst.get(E, X));
	tMapEy.Transfer(src.get(E, Y), dst.get(E, Y));
	tMapEz.Transfer(src.get(E, Z), dst.get(E, Z));
	tMapHx.Transfer(src.get(H, X), dst.get(H, X));
	tMapHy.Transfer(src.get(H, Y), dst.get(H, Y));
	tMapHz.Transfer(src.get(H, Z), dst.get(H, Z));
}

NearToFarFieldReqs::NearToFarFieldReqs(
	const NearToFarFieldProbe& p, const DG_FECollection* fec, FiniteElementSpace& fes, Fields& global) :
	ntff_smsh_{ NearToFarFieldSubMesher(*fes.GetMesh(), fes, buildSurfaceMarker(p.tags, fes)) },
	sfes_{ std::make_unique<FiniteElementSpace>(ntff_smsh_.getSubMesh(), fec) },
	fields_{ Fields(*sfes_) },
	gFields_{ global },
	tMaps_{ TransferMaps(gFields_, fields_) }
{
	updateFields();
}

}
