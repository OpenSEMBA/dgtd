#include "ProbesManager.h"

namespace maxwell {

using namespace mfem;

ParaViewDataCollection ProbesManager::buildParaviewDataCollectionInfo(const ExporterProbe& p, Fields& fields) const
{
	ParaViewDataCollection pd{ p.name, fes_.GetMesh()};
	pd.SetPrefixPath("ParaView");
	
	pd.RegisterField("Ex", &fields(E,X));
	pd.RegisterField("Ey", &fields(E,Y));
	pd.RegisterField("Ez", &fields(E,Z));
	pd.RegisterField("Hx", &fields(H,X));
	pd.RegisterField("Hy", &fields(H,Y));
	pd.RegisterField("Hz", &fields(H,Z));
	
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

	for (auto& p : probes.exporterProbes) {
		p.t_final = opts.t_final;
	}
	
	for (const auto& p : probes.pointProbes) {
		pointProbesCollection_.emplace(&p, buildPointProbeCollectionInfo(p, fields));
	}
}

const PointProbe& ProbesManager::getPointProbe(const std::size_t i) const
{
	assert(i < probes.pointProbes.size());
	return probes.pointProbes[i];
}

const GridFunction& getFieldView(const PointProbe& p, const mfem::FiniteElementSpace& fes, Fields& fields)
{
	switch (p.getFieldType()) {
	case FieldType::E:
		return fields(E, p.getDirection());
	case FieldType::H:
		return fields(H, p.getDirection());
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

void ProbesManager::updateProbe(ExporterProbe& p, double time)
{
	if (cycle_ % p.visSteps != 0 && p.t_final != time) {
		return;
	}

	auto it{ exporterProbesCollection_.find(&p) };
	assert(it != exporterProbesCollection_.end());
	auto& pd{ it->second };

	pd.SetCycle(cycle_);
	pd.SetTime(time);
	pd.Save();
}

void ProbesManager::updateProbe(PointProbe& p, double time)
{
	const auto& it{ pointProbesCollection_.find(&p) };
	assert(it != pointProbesCollection_.end());
	const auto& pC{ it->second };
	p.addFieldToMovie(
		time, 
		pC.field.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP)
	);
}

void ProbesManager::updateProbes(double time)
{
	for (auto& p: probes.exporterProbes) {
		updateProbe(p, time);
	}
	
	for (auto& p : probes.pointProbes) {
		updateProbe(p, time);
	}

	cycle_++;
}
}