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

ProbesManager::ProbesManager(Probes pIn, const mfem::FiniteElementSpace& fes, Fields& fields, const SolverOptions& opts) :
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

	finalTime_ = opts.finalTime;
}

void ProbesManager::initNearToFarFieldProbeDataCollection(NearToFarFieldProbe& p, Fields& fields)
{
		nearToFarFieldProbesCollection_.emplace(&p, buildNearToFarFieldDataCollectionInfo(p, fields));
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

DataCollection ProbesManager::buildNearToFarFieldDataCollectionInfo(const NearToFarFieldProbe& p, Fields& fields)
{
	DataCollection res{ p.name };
	res.SetPrefixPath(p.name);
	res.RegisterField("Ex", &fields.get(E, X));
	res.RegisterField("Ey", &fields.get(E, Y));
	res.RegisterField("Ez", &fields.get(E, Z));
	res.RegisterField("Hx", &fields.get(H, X));
	res.RegisterField("Hy", &fields.get(H, Y));
	res.RegisterField("Hz", &fields.get(H, Z));

	return res;
}

void ProbesManager::updateProbe(ExporterProbe& p, Time time)
{
	if (abs(time - finalTime_) >= 1e-3){
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

void updateFieldsTargets(DataCollection& dc, Fields& fields)
{
	dc.DeregisterField("Ex");
	dc.RegisterField("Ex", &fields.get(E, X));
	dc.DeregisterField("Ey");
	dc.RegisterField("Ey", &fields.get(E, Y));
	dc.DeregisterField("Ez");
	dc.RegisterField("Ez", &fields.get(E, Z));
	dc.DeregisterField("Hx");
	dc.RegisterField("Hx", &fields.get(H, X));
	dc.DeregisterField("Hy");
	dc.RegisterField("Hy", &fields.get(H, Y));
	dc.DeregisterField("Hz");
	dc.RegisterField("Hz", &fields.get(H, Z));
}

FiniteElementSpace reassembleFES(NearToFarFieldProbe& p)
{
	std::string meshDir(p.name + "/" + p.name);
	auto m{ Mesh::LoadFromFile(meshDir) };

	std::string fesDir(p.name + "/" + "fes.txt");
	std::filebuf fb;
	fb.open(fesDir, std::ios::in);
	std::istream iStr(&fb);
	FiniteElementSpace* tfes{};
	auto fec{ tfes->Load(&m, iStr) };

	auto fes{ FiniteElementSpace(&m,fec) };

	return fes;
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

void ProbesManager::updateNearToFarFieldProbe(NearToFarFieldProbe& p, Time time, Fields& fields)
{
	if (abs(time - finalTime_) >= 1e-3) {
		if (cycle_ % p.steps != 0) {
			return;
		}
	}

	auto it{ nearToFarFieldProbesCollection_.find(&p) };
	assert(it != nearToFarFieldProbesCollection_.end());
	auto& dc{ it->second };


	auto fes{ reassembleFES(p) };
	updateFieldsTargets(dc, buildFieldsForProbe(fields, fes));

	dc.SetCycle(cycle_);
	dc.SetTime(time);
	dc.Save();
}

void ProbesManager::updateProbes(Time t, Fields& fields)
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
		updateNearToFarFieldProbe(p, t, fields);
	}

	cycle_++;
}

}