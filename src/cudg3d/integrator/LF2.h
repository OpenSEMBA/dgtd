#pragma once

#include "TimeIntegrator.h"

namespace SEMBA::cudg3d::integrator {

class LF2 : public TimeIntegrator {
public:
	LF2(const dg::Evolution&, const Options&);
//	void setSolver(DG* solver);
//	void timeIntegrate(const Math::Real time) const;
//protected:
//	size_t getNumOfIterationsPerBigTimeStep(const size_t e) const;
//private:
//	static const size_t nStages = 2;
//	size_t getNStages() const;
//	Math::Real getMaxTimeRatio() const;
//	void LTSupdateFieldsElectric(
//	  const Math::Real localTime,
//	  const Math::Real localdt,
//	  const size_t tier) const;
//	void LTSupdateFieldsMagnetic(
//	  const Math::Real localTime,
//	  const Math::Real localdt,
//	  const size_t tier) const;
//	void updateFields(
//	  const size_t e1,
//	  const size_t e2,
//	  const Math::Real localTime,
//	  const Math::Real rkdt) const;
//	void addRHSToFieldsElectric(
//	  const size_t e1,
//	  const size_t e2,
//	  const Math::Real rkdt) const;
//	void addRHSToFieldsMagnetic(
//	  const size_t e1,
//	  const size_t e2,
//	  const Math::Real rkdt) const;
};

}