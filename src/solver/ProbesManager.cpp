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
	pd.SetLevelsOfDetail(order);
	order > 0 ? pd.SetHighOrderOutput(true) : pd.SetHighOrderOutput(false);
	
	pd.SetDataFormat(VTKFormat::BINARY);

	return pd;
}

ProbesManager::ProbesManager(Probes pIn, const mfem::FiniteElementSpace& fes, Fields& fields, const SolverOptions& opts) :
	probes{pIn},
	fes_{fes}
{
	for (const auto& p: probes.exporterProbes) {
		exporterProbesCollection_.emplace(&p, buildParaviewDataCollectionInfo(p, fields));
	}

	for (const auto& p : probes.pointProbes) {
		pointProbesCollection_.emplace(&p, buildPointProbeCollectionInfo(p, fields));
	}

	//for (const auto& p: probes.energyProbes) {
	//	energyProbesCollection_.emplace(&p, buildEnergyProbeCollectionInfo(fes_, fields));
	//}
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

//const EnergyProbe& ProbesManager::getEnergyProbe(const std::size_t i) const
//{
//	assert(i < probes.energyProbes.size());
//	return probes.energyProbes[i];
//}

const GridFunction& getFieldView(const PointProbe& p, const FiniteElementSpace& fes, Fields& fields)
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
		getFieldView(p, fes_, fields)
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
		fields
	};
}

//ProbesManager::EnergyProbeCollection
//ProbesManager::buildEnergyProbeCollectionInfo(const mfem::FiniteElementSpace& fes, const Fields& fields) const
//{
//
//	return {
//		fes,
//		fields
//	};
//}

void ProbesManager::updateProbe(ExporterProbe& p, Time time)
{
	if (cycle_ % p.visSteps != 0) {
		return;
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
	for (auto ft : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			p.addFieldToMovies(
				ft,
				d,
				time,
				pC.fields.get(ft,d).GetValue(pC.fesPoint.elementId, pC.fesPoint.iP)
			);
		}
	}	
}

//void ProbesManager::updateProbe(EnergyProbe& p, Time time)
//{
//	if (cycle_ % p.visSteps != 0) {
//		return;
//	}
//
//	auto& it{ energyProbesCollection_.find(&p) };
//	assert(it != energyProbesCollection_.end());
//	auto& pC{ it->second };
//
//	p.addFieldsToMovie(
//		time,
//		pC.fields
//	);
//}

void ProbesManager::updateProbes(Time t)
{
	for (auto& p : probes.exporterProbes) {
		updateProbe(p, t);
	}
	
	for (auto& p : probes.pointProbes) {
		updateProbe(p, t);
	}

	//for (auto& p : probes.energyProbes) {
	//	updateProbe(p, t);
	//}

	cycle_++;
}
}