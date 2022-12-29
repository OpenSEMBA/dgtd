#include "LSERK4.h"

namespace SEMBA::dgtd::integrator {

const std::size_t nStages{ 5 };
const std::array<Math::Real, nStages> rka{
    0.0,
    -567301805773.0/1357537059087.0,
    -2404267990393.0/2016746695238.0,
    -3550918686646.0/2091501179385.0,
    -1275806237668.0/842570457699.0
};

const std::array<Math::Real, nStages> rkb{
    1432997174477.0/9575080441755.0,
    5161836677717.0/13612068292357.0,
    1720146321549.0/2090206949498.0,
    3134564353537.0/4481467310338.0,
    2277821191437.0/14882151754819.0
};

const std::array<Math::Real, nStages> rkc{
    0.0,
    1432997174477.0/9575080441755.0,
    2526269341429.0/6820363962896.0,
    2006345519317.0/3224310063776.0,
    2802321613138.0/2924317926251.0 
};

LSERK4::LSERK4(const dg::Evolution& dg, const Options& opts):
    TimeIntegrator{opts}
{
    //cfl_ = opts.;
    //buildRKConstants();
    //useMaxStageSizeForLTS = arg->isUseMaxStageSizeForLTS();
    //init(mesh, pmGroup, arg);
}

//void LSERK4::timeIntegrate(const Math::Real time) const {
//    assert(solver != NULL);
//    Math::Real dt = getMaxDT();
//    if (doLTS) {
//        LTSTimeIntegration(time, time, dt, getNTiers() - 1);
//    } else  {
//        Math::Real localTime;
//        const size_t nStages = getNStages();
//        size_t nK = solver->nK;
//        for (size_t s = 0; s < nStages; s++) {
//            localTime = time + dt * getRKC(s);
//            updateResiduals(0, nK, getRKA(s), localTime, dt);
//            solver->updateFieldsWithRes(0, nK, getRKB(s));
//        }
//    }
//}
//
//size_t LSERK4::getNumOfIterationsPerBigTimeStep(
//        const size_t e) const {
//    size_t nTiers = getNTiers();
//    size_t nStages = getNStages();
//    size_t tier = timeTierList_(e,1);
//    size_t stage = timeTierList_(e,2);
//    size_t iter = 0;
//    if (tier == 0) {
//        iter = (nTiers - tier) * nStages;
//    } else {
//        iter = (nTiers - tier) * (nStages + nStages - (stage+1));
//    }
//    return iter;
//}
//
//size_t
//LSERK4::getNStages() const {
//    return nStages;
//}
//
//Math::Real
//LSERK4::getRKA(const size_t s) const {
//    return rka[s];
//}
//
//Math::Real
//LSERK4::getRKB(const size_t s) const {
//    return rkb[s];
//}
//
//Math::Real
//LSERK4::getRKC(const size_t s) const {
//    return rkc[s];
//}
//
//Math::Real
//LSERK4::getStageSize(const size_t s) const {
//    return stageSize[s];
//}
//
//void
//LSERK4::buildRKConstants() {
//    for (size_t i = 1; i <= nStages; i++) {
//        stageSize[i-1] = rkc[i] - rkc[i-1];
//    }
//    stageSize[nStages - 1] = (Math::Real) 1.0 - rkc[nStages-1];
//}
//
//Math::Real
//LSERK4::getMaxStageSize() const {
//    DynMatrix<Math::Real> aux(nStages,1);
//    for (size_t i = 0; i < nStages; i++) {
//        aux(i,0) = stageSize[i];
//    }
//    aux.sortRows(0,0);
//    return aux(nStages-1, 0);
//}
//
//Math::Real
//LSERK4::getMaxTimeRatio() const {
//    if (useMaxStageSizeForLTS) {
//        return getMaxStageSize();
//    } else {
//        return ((Math::Real)1.0 / (Math::Real) nStages);
//    }
//}
//
//
//void
//LSERK4::LTSTimeIntegration(
//        const Math::Real time,
//        Math::Real localTime,
//        Math::Real localdt,
//        const size_t tier) const {
//    size_t fK, lK;
//    const size_t nStages = getNStages();
//    for (size_t s = 0; s < nStages; s++) {
//        // Determines range of cells belonging to this tier and stage.
//        if (tier == getNTiers()-1) {
//            fK = getRange(tier, 0).first;
//            lK = getRange(tier, nStages-1).second;
//        } else {
//            fK = getRange(tier, 0).first;
//            if (s == nStages-1) {
//                lK = getRange(tier, s).second;
//            } else {
//                lK = getRange(tier+1, 3-s).second;
//            }
//        }
//        // Updates RHS in current tier.
//        localTime = time + localdt * getRKC(s);
//        Math::Real rkdt = localdt * getStageSize(s);
//        // Recursively calls this function.
//        if (tier > 0) {
//            size_t lSaved = getRange(tier, nStages-2).second;
//            solver->LTSSaveFieldsAndResidues(fK,lSaved);
//            LTSTimeIntegration(time, localTime, rkdt, tier-1);
//            solver->LTSLoadFieldsAndResidues(fK,lSaved);
//        }
//        // Updates field and residue data in current tier.
//        updateResiduals(fK, lK, getRKA(s), localTime, localdt);
//        solver->updateFieldsWithRes(fK,lK, getRKB(s));
//    }
//}
//
//void
//LSERK4::updateResiduals(
//        const size_t e1,
//        const size_t e2,
//        const Math::Real rka,
//        const Math::Real localTime,
//        const Math::Real localdt) const {
//    solver->computeRHS(e1,e2,localTime,localdt);
//    solver->addRHSToRes(e1,e2,rka,localdt);
//}

}