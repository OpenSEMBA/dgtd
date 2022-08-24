#include "ProbesManager.h"

namespace maxwell {

using namespace mfem;

ParaViewDataCollection ProbesManager::buildParaviewDataCollection(const ExporterProbe& p, Fields& fields) const
{
	ParaViewDataCollection pd{ p.name, fes_.GetMesh()};
	pd.SetPrefixPath("ParaView");
	
	pd.RegisterField("Ex", &fields.E[X]);
	pd.RegisterField("Ey", &fields.E[Y]);
	pd.RegisterField("Ez", &fields.E[Z]);
	pd.RegisterField("Hx", &fields.H[X]);
	pd.RegisterField("Hy", &fields.H[Y]);
	pd.RegisterField("Hz", &fields.H[Z]);

	const auto order{ fes_.GetMaxElementOrder() };
	pd.SetLevelsOfDetail(order);
	order > 0 ? pd.SetHighOrderOutput(true) : pd.SetHighOrderOutput(false);
	
	pd.SetDataFormat(VTKFormat::BINARY);

	return pd;
}

ProbesManager::ProbesManager(Probes probes, const mfem::FiniteElementSpace& fes, Fields& fields) :
	probes_{probes},
	fes_{fes}
{
	for (const auto& p: probes_.exporterProbes) {
		exporterProbesCollection_.emplace(&p, buildParaviewDataCollection(p, fields));
	}
	
	for (const auto& p : probes_.pointsProbes) {
		pointProbesCollection_.emplace(&p, buildPointsProbeCollection(p, fields));
	}
}

const PointsProbe& ProbesManager::getPointsProbe(const std::size_t i) const
{
	assert(i < probes_.pointsProbes.size());
	return probes_.pointsProbes[i];
}

const GridFunction& getFieldView(const PointsProbe& p, Fields& fields)
{
	switch (p.getFieldType()) {
	case FieldType::E:
		return fields.E[p.getDirection()];
	case FieldType::H:
		return fields.H[p.getDirection()];
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

ProbesManager::PointsProbeCollection
ProbesManager::buildPointsProbeCollection(const PointsProbe& p, Fields& fields) const
{
	std::vector<FESPoint> fesPoints;
	for (const auto& point : p.getPoints()) {
		Array<int> elemIdArray;
		Array<IntegrationPoint> integPointArray;
		fes_.GetMesh()->FindPoints(toDenseMatrix(point), elemIdArray, integPointArray);
		assert(elemIdArray.Size() == 1);
		assert(integPointArray.Size() == 1);
		fesPoints.push_back({ elemIdArray[0], integPointArray[0] });
	}
	
	return { 
		fesPoints, 
		getFieldView(p, fields)
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

void ProbesManager::updateProbe(PointsProbe& p, double time)
{
	const auto& it{ pointProbesCollection_.find(&p) };
	assert(it != pointProbesCollection_.end());
	const auto& pC{ it->second };

	FieldFrame frame;
	for (const auto& fesPoint: pC.fesPoints) {
		frame.push_back( pC.field.GetValue(fesPoint.elementId, fesPoint.iP) );
	}
	
	p.addFrame(time, frame);
}

void ProbesManager::updateProbes(double time)
{
	if (cycle_ % probes_.visSteps == 0) {
		for (auto& p: probes_.exporterProbes) {
			updateProbe(p, time);
		}
		for (auto& p : probes_.pointsProbes) {
			updateProbe(p, time);
		}
	}
	cycle_++;
}
}