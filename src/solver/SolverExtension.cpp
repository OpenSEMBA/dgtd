#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

namespace maxwell
{

using namespace mfem;

InteriorFaceConnectivityMaps getGlobalNodeID(const InteriorFaceConnectivityMaps& local_dof_ids, const GlobalConnectivity& global)
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

void SBCSolver::findDoFPairs(Model& model, ParFiniteElementSpace& fes)
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

GeomTagToMaterial getSBCSolverGeomTagToMaterialFromGlobal(Model& g_model)
{
    GeomTagToMaterial res;
    const auto& tag2bdr = g_model.getGeomTagToIntBoundaryCond();
    for(const auto& [tag, cond] : tag2bdr){
        if (cond == BdrCond::SBC){
            const auto sbc_material = g_model.getGeomTagToBoundaryMaterial().at(tag);
            res.emplace(1, sbc_material);
        }
    }
    return res;
}

SBCSolver::SBCSolver(Model& g_model, ParFiniteElementSpace& g_fes, const SBCProperties& sbcp) :
sbcp_(sbcp)
{

    auto mesh  = Mesh::MakeCartesian1D(sbcp_.num_of_segments + 2, sbcp_.material_width + 2 * (sbcp_.material_width / double(sbcp_.num_of_segments)));
    auto pmesh = ParMesh(MPI_COMM_WORLD, mesh);
    auto fec   = DG_FECollection(sbcp_.order, 1, BasisType::GaussLobatto);
    auto pfes  = ParFiniteElementSpace(&pmesh, &fec);

    Model model(mesh, GeomTagToMaterialInfo(getSBCSolverGeomTagToMaterialFromGlobal(g_model), GeomTagToBoundaryMaterial{}), GeomTagToBoundaryInfo{});
    Probes pr;
    Sources src;
    EvolutionOptions eopts;
    eopts.order = sbcp_.order;
    ProblemDescription pd(model, pr, src, eopts);
    DGOperatorFactory<ParFiniteElementSpace> dgops(pd, pfes);
    auto global_operator = dgops.buildGlobalOperator();
    global_operator.get()->Print(std::cout);

    findDoFPairs(g_model, g_fes);
    
    // assignEvolutionOperator();
    // model_ = Model(*mesh_, GeomTagToMaterialInfo(getSBCSolverGeomTagToMaterialFromGlobal(g_model), GeomTagToBoundaryMaterial{}));


    
}

void SBCSolver::assignGlobalFields(const Fields<ParFiniteElementSpace,ParGridFunction>* g_fields)
{
    global_fields_ = g_fields;
}

// void SBCSolver::loadFieldValues(const FieldType f, const Direction d, const NbrPairs& vals)
// {
//     const auto& ghost_interval = sbcp_.order + 1;
//     for(auto v = 0; v < ghost_interval; v++){
//         this->sbc_fields_.get(f, d)[v] = vals.first;
//         this->sbc_fields_.get(f, d)[sbc_fields_.get(f, d).Size() - v] = vals.second;
//     }
// }

// NbrPairs SBCSolver::getFieldValues(const FieldType f, const Direction d)
// {
//     const auto& ghost_interval = sbcp_.order + 1;
//     return {this->sbc_fields_.get(f,d)[sbcp_.order + 2], this->sbc_fields_.get(f,d)[sbc_fields_.get(f,d).Size() - (sbcp_.order + 2)]};
// }

SBCTimeDependentOperator::SBCTimeDependentOperator(Model& model, ParFiniteElementSpace& fes) :
model_(model),
fes_(fes)
{
    Probes pr;
    Sources src;
    EvolutionOptions eopts;
    ProblemDescription pd(model_, pr, src, eopts);
    DGOperatorFactory<ParFiniteElementSpace> dgops(pd, fes_);

    sbc_operator_ = dgops.buildGlobalOperator();
    
}

}