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

#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>			// Stream I/O.
#include <cmath>
#include <vector>
#include <utility>

using namespace std;

#include "SmbData.h"
#include "CellGroup.h"

class BoundaryCondition {
public:
    BoundaryCondition();
    BoundaryCondition(const CellTet<ORDER_N>* cell, UInt face);
    virtual ~BoundaryCondition();
    bool hasSameBoundary(const BoundaryCondition& other) const;
    virtual BoundaryCondition& operator=(const BoundaryCondition& rhs);
    virtual void printInfo() const;
    const CellTet<ORDER_N>* getCell() const;
    UInt getFace() const;
    Face getCellFace() const;
private:
    const CellTet<ORDER_N>* cell_;
    UInt face_;
};

class EMSourceBC : public BoundaryCondition {
public:
    EMSourceBC();
    virtual ~EMSourceBC();
    EMSourceBC(const CellTet<ORDER_N>* e, const UInt f, const EMSourceBase* bc);
    EMSourceBC& operator=(const EMSourceBC& rhs);
    void check() const;
    const EMSourceBase* getCondition() const;

private:
    const EMSourceBase* em_;
};

class PhysicalModelBC : public BoundaryCondition {
public:
    //
    PhysicalModelBC();
    PhysicalModelBC(
            const CellTet<ORDER_N>*,
            UInt face,
            const PMPredefined* bc);
    virtual ~PhysicalModelBC();
    PhysicalModelBC& operator=(const PhysicalModelBC& rhs);
    const PMPredefined* getCondition() const;
private:
    const PMPredefined* pm_;
};

class SurfaceImpedanceBC : public BoundaryCondition {
public:
    SurfaceImpedanceBC();
    virtual ~SurfaceImpedanceBC();
    SurfaceImpedanceBC(
            const CellTet<ORDER_N>* cell,
            const UInt face,
            const CellTet<ORDER_N>* cellD,
            const UInt faceD,
            const PMSurfaceSIBC* cond);
    SurfaceImpedanceBC& operator=(const SurfaceImpedanceBC &rhs);
    const CellTet<ORDER_N>* getCellD() const;
    UInt getFaceD() const;
    const PMSurfaceSIBC* getCondition() const;

private:
    const CellTet<ORDER_N>* cellD_;
    UInt faceD_;
    const PMSurfaceSIBC* sibc_;
};


#endif /* BOUNDARYCONDITION_H_ */
