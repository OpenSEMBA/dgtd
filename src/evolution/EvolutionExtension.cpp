#include "EvolutionExtension.h"

namespace maxwell
{

using namespace mfem;

static InteriorFaceConnectivityMaps getGlobalNodeID(const InteriorFaceConnectivityMaps& local_dof_ids, const GlobalConnectivity& global)
	{
		InteriorFaceConnectivityMaps res;
		res.first.resize(local_dof_ids.first.size());
		res.second.resize(local_dof_ids.second.size());
		for (auto v{ 0 }; v < res.first.size(); v++) {
			res.first[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(local_dof_ids.first[v], local_dof_ids.second[v])));
			res.second[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(local_dof_ids.second[v], local_dof_ids.first[v])));
		}
		return res;
	}

void SBCManager::findDoFPairs(Model& model, FiniteElementSpace& fes)
{
    auto attMap{ mapOriginalAttributes(*fes.GetMesh()) };
    auto fec = dynamic_cast<const DG_FECollection*>(fes.FEColl());
    GlobalConnectivity global = assembleGlobalConnectivityMap(*fes.GetMesh(), fec);
    auto sbc_marker = model.getMarker(BdrCond::SBC, true);
    for (auto b = 0; b < model.getMesh().GetNBE(); b++){
        if (sbc_marker[model.getConstMesh().GetBdrAttribute(b) - 1] == 1) {
            const FaceElementTransformations* faceTrans;
            fes.GetMesh()->FaceIsInterior(fes.GetMesh()->GetFaceElementTransformations(fes.GetMesh()->GetBdrElementFaceIndex(b))->ElementNo) ? faceTrans = fes.GetMesh()->GetInternalBdrFaceTransformations(b) : faceTrans = fes.GetMesh()->GetBdrFaceTransformations(b);
            auto twoElemSubMesh{ assembleInteriorFaceSubMesh(*fes.GetMesh(), *faceTrans, attMap) };
            FiniteElementSpace subFES(&twoElemSubMesh, fec);
            auto node_pair_global{ getGlobalNodeID(buildConnectivityForInteriorBdrFace(*faceTrans, fes, subFES), global)};
            for (auto p = 0; p < node_pair_global.first.size(); p++){
                dof_pairs_.emplace_back(node_pair_global.first[p], node_pair_global.second[p]);
            }
        }
    }
}

SBCManager::SBCManager(Model& model, FiniteElementSpace& fes, const SBCProperties& sbcp)
{
    findDoFPairs(model, fes);

    mesh_ = std::make_unique<Mesh>(Mesh::MakeCartesian1D(sbcp.num_of_segments, sbcp.material_width));
    fec_ = std::make_unique<DG_FECollection>(sbcp.order, 1, BasisType::GaussLobatto);
    fes_ = std::make_unique<FiniteElementSpace>(mesh_.get(), fec_.get());
}

}