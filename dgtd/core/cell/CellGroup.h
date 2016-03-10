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
 * CellGroup.h
 *
 *  Created on: Aug 29, 2012
 *      Author: luis
 */

#ifndef CELLGROUP_H_
#define CELLGROUP_H_

#include "../../dgtd/core/CellTet.h"
#include "../../dgtd/core/Ordering.h"

class CellGroup : public Ordering {
public:
	vector<CellTet<ORDER_N>*> cell;
	vector<CellTet4<ORDER_N> > linTet;
	vector<CellTet10<ORDER_N> > quadTet;
	CellGroup(const Mesh::Volume& mesh, const PMGroup& pMGroup);
	~CellGroup();
	const CellTet<ORDER_N>* operator()(const UInt i) const;
	const CellTet<ORDER_N>* getPtrToCell(const VolR* elem) const;
	const CellTet<ORDER_N>* getPtrToCellWithId(const ElemId&) const;
private:
	UInt cellOffsetId;
	void buildNodalMaps(const Connectivities& map);
	void check(const Connectivities& map) const;
	void checkNodalMaps(const Connectivities& map) const;
};

#endif /* CELLGROUP_H_ */
