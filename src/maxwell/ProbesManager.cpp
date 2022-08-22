#include "ProbesManager.h"

namespace maxwell {

using namespace mfem;

ParaViewDataCollection buildParaviewDataCollection(
	FiniteElementSpace* fes, FieldViews& fields)
{
	ParaViewDataCollection pd{ "MaxwellView", fes->GetMesh()};
	pd.SetPrefixPath("ParaView");
	
	pd.RegisterField("Ex", &(*fields.E)[X]);
	pd.RegisterField("Ey", &(*fields.E)[Y]);
	pd.RegisterField("Ez", &(*fields.E)[Z]);
	pd.RegisterField("Hx", &(*fields.H)[X]);
	pd.RegisterField("Hy", &(*fields.H)[Y]);
	pd.RegisterField("Hz", &(*fields.H)[Z]);

	const auto order{ fes->GetMaxElementOrder() };
	pd.SetLevelsOfDetail(order);
	order > 0 ? pd.SetHighOrderOutput(true) : pd.SetHighOrderOutput(false);
	
	pd.SetDataFormat(VTKFormat::BINARY);

	return pd;
}

ProbesManager::ProbesManager(Probes probes, const mfem::FiniteElementSpace* fes, const FieldViews& fields) :
	probes_{probes}
{
	for (const auto& p: probes_.exporterProbes) {
		buildParaviewDataCollection(fes, fields);
	}
	
	for (const auto& p : probes_.pointsProbes) {
		probesToFES_.emplace(&p, buildElemAndIntegrationPointArrays(p));
	}
}

const PointsProbe* ProbesManager::getPointsProbe(const std::size_t i) const
{
	assert(i < probes_.pointsProbes.size());
	return &probes_.pointsProbes[i];
}

ProbesManager::PointsProbeInFES
buildElemAndIntegrationPointArrays(const PointsProbe& p, const mfem::FiniteElementSpace* fes)
{
	DenseMatrix physPoints{ p.getIntegPointMat() };
	Array<int> elemIdArray;
	Array<IntegrationPoint> integPointArray;
	std::pair<Array<int>, Array<IntegrationPoint>> res;
	fes->GetMesh()->FindPoints(physPoints, elemIdArray, integPointArray);
	return { elemIdArray, buildIntegrationPointsSet(integPointArray) };
}

const std::vector<std::vector<IntegrationPoint>> buildIntegrationPointsSet(const Array<IntegrationPoint>& ipArray)
{
	std::vector<IntegrationPoint> aux;
	aux.resize(model_.getConstMesh().Dimension());
	IntegrationPointsSet res;
	res.resize(ipArray.Size(), aux);
	for (int i = 0; i < ipArray.Size(); i++) {
		switch (fes_.GetMesh()->Dimension()) {
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

FieldFrame ProbesManager::getFieldForPointsProbe(const PointsProbe& p) const
{
	const auto& it{ probeIntegrationPoints_.find(&p) };
	assert(it != probeIntegrationPoints_.end());
	const auto& pPIP{ it->second };

	FieldFrame res;
	res.resize(pPIP.elemIds.Size());
	for (int j = 0; j < pPIP.elemIds.Size(); j++) {
		assert(pPIP.elemIds.Size() == pPIP.integPointSet.size());
		for (int d = X; d <= Z; d++) {
			const auto& dir{ p.getDirection() };
			switch (p.getFieldType()) {
			case FieldType::E:
				res[j][dir] = E_[dir].GetValue(pPIP.elemIds[j], pPIP.integPointSet[j][d]);
				break;
			case FieldType::H:
				res[j][dir] = H_[dir].GetValue(pPIP.elemIds[j], pPIP.integPointSet[j][d]);
				break;
			}
		}
	}
	return res;
}


void ProbesManager::storeInitialVisualizationValues()
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

void ProbesManager::updateProbes(double time, bool done, int cycle)
{
	storeInitialVisualizationValues();
	for (auto& p : probes_.pointsProbes) {
		p.addFrame(time, getFieldForPointsProbe(p));
	}
	if (done || cycle % probes_.vis_steps == 0) {
		for (auto& p : probes_.getPointsProbes()) {
			p.addFrame(time, getFieldForPointsProbe(p));
		}
		for (int i = 0; i < probes_.getExporterProbes().size(); i++) {
			if (probes_.getExporterProbes().at(i).type == ExporterProbe::Type::Paraview) {
				pd_->SetCycle(cycle);
				pd_->SetTime(time);
				pd_->Save();
				break;
			}
		}
	}

}
}