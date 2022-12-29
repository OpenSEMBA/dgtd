#pragma once

#include "TimeIntegrator.h"

namespace SEMBA::cudg3d::integrator {

class LSERK4 : public TimeIntegrator {
public:
    LSERK4(const dg::Evolution&, const Options&);
//    void timeIntegrate(const Math::Real time) const;
//protected:
//    size_t getNumOfIterationsPerBigTimeStep(
//            const size_t e) const;
//private:
//    bool useMaxStageSizeForLTS;
//    static const size_t nStages = 5;
//    static const Math::Real rka[nStages];
//    static const Math::Real rkb[nStages];
//    static const Math::Real rkc[nStages];
//    Math::Real stageSize[nStages];
//    Math::Real getMaxStageSize() const;
//    Math::Real getMaxTimeRatio() const;
//    void buildRKConstants();
//    void LTSTimeIntegration(
//            const Math::Real time,
//            const Math::Real localTime,
//            const Math::Real localdt,
//            const size_t tier) const;
//    void updateResiduals(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real rka,
//            const Math::Real localTime,
//            const Math::Real localdt) const;
//    size_t getNStages() const;
//    Math::Real getRKA(const size_t s) const;
//    Math::Real getRKB(const size_t s) const;
//    Math::Real getRKC(const size_t s) const;
//    Math::Real getStageSize(const size_t s) const;
};

}