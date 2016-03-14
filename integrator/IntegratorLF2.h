// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
/*
 * SolverLeapfrog.h
 *
 *  Created on: Feb 22, 2013
 *      Author: luis
 */

#ifndef INTEGRATORLF2_H_
#define INTEGRATORLF2_H_

#include "../../dgtd/integrator/Integrator.h"

#ifdef LEAPFROG_ORDER
#ifndef SOLVER_IGNORE_DISPERSIVES
	#error Dispersive materials have not been implemented for LF.
#endif
#endif

class IntegratorLF2 : public Integrator {
public:
	IntegratorLF2();
	virtual ~IntegratorLF2();
	IntegratorLF2(
	 const Mesh::Volume& mesh,
	 const PMGroup& pmGroup,
	 const OptionsSolverDGTD* arg);
	void
	 setSolver(DG* solver);
	void
	 timeIntegrate(const Real time) const;
protected:
	size_t
 	 getNumOfIterationsPerBigTimeStep(
      const size_t e) const;
private:
	static const size_t nStages = 2;
	size_t
	 getNStages() const;
	Real
	 getMaxTimeRatio() const;
	void
	 LTSupdateFieldsElectric(
	  const Real localTime,
	  const Real localdt,
	  const size_t tier) const;
	void
	 LTSupdateFieldsMagnetic(
	  const Real localTime,
	  const Real localdt,
	  const size_t tier) const;
	void
	 updateFields(
	  const size_t e1,
	  const size_t e2,
	  const Real localTime,
	  const Real rkdt) const;
	void
	 addRHSToFieldsElectric(
	  const size_t e1,
	  const size_t e2,
	  const Real rkdt) const;
	void
	 addRHSToFieldsMagnetic(
	  const size_t e1,
	  const size_t e2,
	  const Real rkdt) const;
};

#endif /* SOLVERLEAPFROG_H_ */
