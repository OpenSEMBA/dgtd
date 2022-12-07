#pragma once

#include "Source.h"

namespace SEMBA::dgtd::dg::source {

class Waveport : public Source {
public:
    //	Math::Real
    //     getNumericalGammaMGauss(
    //      const Math::Real time,
    //      const Math::Real minDT,
    //      const Math::Real amplitude,
    //      const Math::Real delay,
    //      const Math::Real spread,
    //      const Math::Real kcm) const;
    bool checkNormalsAreEqual(
            const vector<pair<size_t,size_t> >& elemFace) const;
protected:
    CVecR3* posTF   = nullptr;
    CVecR3* posTFNB = nullptr;
    CVecR3* posSF   = nullptr;
private:
    Math::Real* gauss = nullptr;
    Math::Real hm = nullptr;
    Math::Real getHm(
            const Math::Real time,
            const Math::Real kcm) const;
};

}
