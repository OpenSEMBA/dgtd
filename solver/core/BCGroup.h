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
 * BCGroup.h
 *
 *  Created on: Jul 8, 2013
 *      Author: luis
 */

#ifndef BCGROUP_H_
#define BCGROUP_H_

#include "BoundaryCondition.h"

using namespace std;

class GroupBoundaryConditions {
public:
    vector<EMSourceBC> embc;
    vector<PhysicalModelBC> pmbc;
    vector<SurfaceImpedanceBC> sibc;
    GroupBoundaryConditions(
            const MeshVolume& mesh,
            const EMSourceGroup& em,
            const PMGroup& pm,
            const CellGroup& cells,
            const Connectivities& map);
    GroupBoundaryConditions& operator=(const GroupBoundaryConditions &rhs);
    vector<const BoundaryCondition*> getPtrsToPEC() const;
    vector<const BoundaryCondition*> getPtrsToPMC() const;
    vector<const BoundaryCondition*> getPtrsToSMA() const;
    vector<const BoundaryCondition*> getPtrsToSIBC() const;
    vector<const BoundaryCondition*> getPtrsToEMSourceBC() const;
    vector<const BoundaryCondition*> getPtrsToBC(const EMSourceBase* pw) const;
    vector<const BoundaryCondition*> getPtrsToBCWithMatId(const MatId id) const;
    void printInfo() const;
private:
    void buildEMSourceBC(
            const MeshVolume& mesh,
            const EMSourceGroup& em,
            const CellGroup& cells);
    void buildPhysicalModelBC(
            const MeshVolume& mesh,
            const PMGroup& pm,
            const CellGroup& cells,
            const Connectivities& map);
    void removeOverlapped();
    vector<BoundaryCondition*> removeCommons(
            const vector<BoundaryCondition*>& low,
            const vector<BoundaryCondition*>& high) const;
    void check() const;
    bool checkOverlapping() const;
    void checkEMSourcesAreSetInVacuum() const;
};

typedef GroupBoundaryConditions BCGroup;

#endif /* BCGROUP_H_ */
