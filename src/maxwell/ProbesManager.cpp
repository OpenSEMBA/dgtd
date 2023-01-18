#include "ProbesManager.h"

namespace maxwell {

using namespace mfem;

ParaViewDataCollection ProbesManager::buildParaviewDataCollectionInfo(const ExporterProbe& p, Fields& fields) const
{
	ParaViewDataCollection pd{ p.name, fes_.GetMesh()};
	pd.SetPrefixPath("ParaView");
	
	switch (fes_.GetMesh()->Dimension()) {
	case 1:
		pd.RegisterField("E", &fields.E1D);
		pd.RegisterField("H", &fields.H1D);
		break;
	default:
		pd.RegisterField("Ex", &fields.E[X]);
		pd.RegisterField("Ey", &fields.E[Y]);
		pd.RegisterField("Ez", &fields.E[Z]);
		pd.RegisterField("Hx", &fields.H[X]);
		pd.RegisterField("Hy", &fields.H[Y]);
		pd.RegisterField("Hz", &fields.H[Z]);
		break;
	}

	const auto order{ fes_.GetMaxElementOrder() };
	pd.SetLevelsOfDetail(order);
	order > 0 ? pd.SetHighOrderOutput(true) : pd.SetHighOrderOutput(false);
	
	pd.SetDataFormat(VTKFormat::BINARY);

	return pd;
}

ProbesManager::ProbesManager(Probes pIn, const mfem::FiniteElementSpace& fes, Fields& fields) :
	probes{pIn},
	fes_{fes}
{
	for (const auto& p: probes.exporterProbes) {
		exporterProbesCollection_.emplace(&p, buildParaviewDataCollectionInfo(p, fields));
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
		switch (fes.GetMesh()->Dimension()) {
		case 1:
			return fields.E1D;
			break;
		default:
			return fields.E[p.getDirection()];
			break;
		}
	case FieldType::H:
		switch (fes.GetMesh()->Dimension()) {
		case 1:
			return fields.H1D;
			break;
		default:
			return fields.H[p.getDirection()];
			break;
		}
	default:
		throw std::runtime_error("Invalid field type.");
	}
}

DenseMatrix toDenseMatrix(const Point& p)
{
	DenseMatrix r{ 1, (int)p.size() };
	for (auto i{ 0 }; i < p.size(); ++i) {
		r(0, i) = p[i];
	}
	return r;
}

ProbesManager::PointProbeCollection
ProbesManager::buildPointProbeCollectionInfo(const PointProbe& p, Fields& fields) const
{
	
	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	fes_.GetMesh()->FindPoints(toDenseMatrix(p.getPoint()), elemIdArray, integPointArray);
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
	if (cycle_ % probes.visSteps == 0) {
		for (auto& p: probes.exporterProbes) {
			updateProbe(p, time);
		}
		for (auto& p : probes.pointProbes) {
			updateProbe(p, time);
		}
	}
	cycle_++;
}
}