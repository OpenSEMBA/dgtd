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

void SBCSolver::findDoFPairs(Model& model, FiniteElementSpace& fes)
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

void SBCSolver::estimateTimeStep()
{
    dt_ = getMinimumInterNodeDistance(*fes_.get()) / std::pow(double(fes_.get()->FEColl()->GetOrder()), 1.5) / physicalConstants::speedOfLight;
}

SBCSolver::SBCSolver(Model& model, FiniteElementSpace& full_model_fes, const SBCProperties& sbcp) :
sbcp_(sbcp),
mesh_(std::make_unique<Mesh>(Mesh::MakeCartesian1D(sbcp_.num_of_segments, sbcp_.material_width))),
fec_(std::make_unique<DG_FECollection>(sbcp_.order, 1, BasisType::GaussLobatto)),
fes_(std::make_unique<FiniteElementSpace>(mesh_.get(), fec_.get())),
sbc_fields_(Fields<FiniteElementSpace,GridFunction>(*fes_.get()))
{
    assignODESolver();
    assignEvolutionOperator();

    findDoFPairs(model, full_model_fes);

    
}

void SBCSolver::assignGlobalFields(const Fields<ParFiniteElementSpace,ParGridFunction>* g_fields)
{
    global_fields_ = g_fields;
}

void SBCSolver::assignODESolver()
{
    switch(sbcp_.implicit_ode){
        case true:
            odeSolver_ = std::make_unique<ImplicitMidpointSolver>();
            break;
        case false:
            odeSolver_ = std::make_unique<RK4Solver>();
            break;
    }
}

void SBCSolver::resetFields()
{
    this->sbc_fields_.allDOFs() = 0.0;
}

std::pair<double, double> SBCSolver::getFieldPairAfterCalculation(const FieldType f, const Direction d)
{
    return {this->sbc_fields_.get(f,d)[0], this->sbc_fields_.get(f,d)[sbc_fields_.get(f,d).Size() - 1]};
}

void SBCSolver::assignEvolutionOperator()
{
    // evolTDO_ = std::make_unique<SBC_TDO>();  // WIP
}

SBCTimeDependentOperator::SBCTimeDependentOperator(Model& model, FiniteElementSpace& fes) :
model_(model),
fes_(fes)
{
    Probes pr;
    Sources src;
    EvolutionOptions eopts;
    ProblemDescription pd(model_, pr, src, eopts);
    DGOperatorFactory<FiniteElementSpace> dgops(pd, fes_);

    sbc_operator_ = dgops.buildGlobalOperator();
    
}

}