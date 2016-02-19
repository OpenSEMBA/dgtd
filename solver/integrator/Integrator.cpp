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
 * Integrator.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: luis
 */

#include "Integrator.h"

Integrator::Integrator() {
    mindt = 0.0;
    nTiers_ = 0;
    tierRange_ = NULL;
    doLTS = true;
    growSmallerTiers = 0;
    maxNumOfTiers = 0;
    cfl_ = 0.0;
    solver = NULL;
}

Integrator::~Integrator() {
    for (UInt i = 0; i < nTiers_; i++) {
        delete tierRange_[i];
    }
}

void Integrator::setSolver(DG* solver_) {
    solver = solver_;
}

UInt Integrator::getNTiers() const {
    return nTiers_;
}

Real Integrator::getMaxDT() const {
    if (doLTS) {
        return (mindt / pow(getMaxTimeRatio(), Real(nTiers_-1)));
    } else {
        return mindt;
    }
}

Real Integrator::getMinDT() const {
    return mindt;
}

vector<vector<ElementId>> Integrator::getTiersIds() const {
    vector<vector<ElementId>> res;
    for (UInt tier = 0; tier < nTiers_; tier++) {
        vector<ElementId> tierIds = getIdsOfTier(tier);
        res.push_back(tierIds);
    }
    return res;
}

vector<vector<ElementId>> Integrator::getStagesIds() const {
    vector<vector<ElementId> > res;
    for (UInt stage = 0; stage < getNStages(); stage++) {
        vector<ElementId> stageIds = getIdsOfStage(stage);
        res.push_back(stageIds);
    }
    return res;
}

vector<vector<ElementId>> Integrator::getPartitionsIds() const {
    //assert(partIds.size() != 0);
    return partIds_;
}

vector<pair<ElementId,Int>> Integrator::getComputationalWeights(
        const MeshVolume* msh) const {
    const Int curlFlops = 1;
    const Int fluxFlops = 0;
    Int flops = curlFlops + fluxFlops;
    GroupElements<const VolR> physVol = msh->elems();
    physVol.removeMatId(MatId(0));
    const UInt nK = physVol.sizeOf<VolR>();
    vector<pair<ElementId,Int>> idWgt;
    idWgt.reserve(nK);
    for (UInt e = 0; e < nK; e++) {
        pair<ElementId,Int> aux;
        ElementId id = getIdOfGlobalRelPos(e);
        aux.first = id;
        aux.second = getNumOfIterationsPerBigTimeStep(e) * flops;
        idWgt.push_back(aux);
    }
    return idWgt;
}

Interval Integrator::getRange(const UInt tier, const UInt stage) const {
    assert(tier < nTiers_);
    assert(stage < getNStages());
    assert(tierRange_ != NULL);
    return tierRange_[tier][stage];
}

UInt Integrator::getNPartitions() const {
    return partIds_.size();
}

void Integrator::partitionate(
        const MeshVolume* msh,
        Comm* comm) {
    cout << " - Getting computational weights... " << flush;
    vector<pair<ElementId,Int> > idWgt = getComputationalWeights(msh);
    cout << "OK" << endl;
    cout << " - Obtaining partition ids... " << flush;
    vector<vector <ElementId>> partId =
            msh->getPartitionsIds(comm->getNumberOfTasks(), idWgt);
    cout << "OK" << endl;
    cout << " - Setting partition sizes... " << flush;
    comm->setPartitionSizes(partId);
    cout << "OK" << endl;
    cout << " - Reordering partitions... " << flush;
    reorder(partId, comm->getLocalOffset(), comm->getLocalSize());
    cout << "OK" << endl;
}

void Integrator::printInfo() const {
    cout << "--- SolverInfo ---" << endl;
    cout << "Min. time step: " << mindt*1E12 << " [ps]" << endl;
    cout << "Max. time step: " << getMaxDT()*1E12 << " [ps]" << endl;
    cout << "Number of tiers: " << nTiers_ << endl;
    if (nTiers_ > 1) {
        for (UInt i = 0; i < nTiers_; i++) {
            cout << "# of Cells in tier " << i << ": "
                    << getNumberOfCellsInTier(i) << endl;
            for (UInt j = 0; j < getNStages(); j++) {
                cout << "--> Range stage #" << j << ": "
                        << getRange(i,j).first << ", "
                        << getRange(i,j).second << endl;
            }
        }
    }
}

void Integrator::init(
        const MeshVolume& mesh,
        const PMGroup& pmGroup,
        const OptionsSolverDGTD* arg) {
    growSmallerTiers = arg->getGrowSmallerTiers();
    maxNumOfTiers = arg->getMaxNumberOfTiers();
    doLTS = arg->isUseLTS();
    if (!doLTS) {
        maxNumOfTiers = 1;
    }
    cout << "- Building tier info... " << flush;
    buildTierInfo(mesh, pmGroup);
    cout << " OK" << endl;
    cout << "- Building Relative Positions of Ids... " << flush;
    buildRelPosOfIds(timeTierList_);
    cout << " OK" << endl;
    cout << "- Building Tier Range... " << flush;
    tierRange_ = buildTierRange(tierRange_, timeTierList_);
    cout <<  "OK" << endl;
}

Real Integrator::getMaxTimeStep(
        const VolR* vol,
        const PhysicalModel* mat) const {
    // Returns the maximum time step allowed for this cell.
    Real fS1 = 0.0;
    for (UInt f = 0; f < vol->numberOfFaces(); f++) {
        Real area = vol->getAreaOfFace(f);
        Real volume = vol->getVolume();
        Real fS2 = area / volume;
        if (fS2 > fS1) {
            fS1 = fS2;
        }
    }
    Real dt = (1.0 / Constants::c0) *  (1.0 / Real(fS1 * ORDER_N * ORDER_N));
    // Checks case of electrical dispersive materials.
    checkMaterialStabilityForDT(mat, dt);
    return (dt * cfl_);
}

UInt Integrator::getNumberOfCellsInTier(const UInt tier) const {
    assert(tier < nTiers_);
    assert(tierRange_ != NULL);
    UInt first = getRange(tier, 0).first;
    UInt last = getRange(tier, getNStages() - 1).second;
    UInt res = last - first;
    return res;
}

vector<ElementId> Integrator::getIdsOfTier(UInt tier) const {
    assert(tier < nTiers_);
    const UInt nK = timeTierList_.nRows();
    vector<ElementId> res;
    res.reserve(nK);
    for (UInt i = 0; i < nK; i++) {
        if (timeTierList_(i,1) == tier) {
            res.push_back(ElementId(timeTierList_(i,0)));
        }
    }
    return res;
}

vector<ElementId> Integrator::getIdsOfStage(UInt stage) const {
    assert(stage < getNStages());
    UInt nK = timeTierList_.nRows();
    vector<ElementId> res;
    res.reserve(nK);
    for (UInt i = 0; i < nK; i++) {
        if (timeTierList_(i,2) == stage) {
            res.push_back(ElementId(timeTierList_(i,0)));
        }
    }
    return res;
}

void Integrator::reorder(
        const vector<vector<ElementId> >& partitionId,
        const UInt localOffset,
        const UInt localSize) {
    partIds_ = partitionId;
    // Reorders timeTierList according to partitions.
    reorderTimeTierList(partitionId);
    // Now aux stores the final ordering considering partitions.
    buildRelPosOfIds(timeTierList_);
    // Builds time tier range
    DynMatrix<UInt> reducedList(localSize,3);
    for (UInt i = 0; i < localSize; i++) {
        for (UInt j = 0; j < 3; j++) {
            reducedList(i,j) = timeTierList_(i + localOffset, j);
        }
    }
    tierRange_ = buildTierRange(tierRange_, reducedList);
}

void Integrator::reorderTimeTierList(
        const vector<vector<ElementId>>& partitionId) {
    UInt nK = timeTierList_.nRows();
    DynMatrix<UInt> aux(nK, 5); // relPos - Ids - Part - Tier - Stage
    for (UInt k = 0; k < nK; k++) {
        aux(k, 0) = k;
        aux(k, 1) = timeTierList_(k, 0);
        aux(k, 2) = 0; // Temporary assignation of partition.
        aux(k, 3) = timeTierList_(k, 1); // Tier assignation.
        aux(k, 4) = timeTierList_(k, 2);
    }
    aux.sortRows_omp(1, 1);
    UInt initId = aux(0, 1);
    UInt nParts = partitionId.size();
    for (UInt p = 0; p < nParts; p++) {
        for (UInt i = 0; i < partitionId[p].size(); i++) {
            UInt id = partitionId[p][i].toUInt();
            aux(id - initId, 2) = p;
        }
    }
    aux.sortRows_omp(2, 4);
    for (UInt i = 0; i < nK; i++) {
        timeTierList_(i, 0) = aux(i, 1);
        timeTierList_(i, 1) = aux(i, 3);
        timeTierList_(i, 2) = aux(i, 4);
    }
}

void Integrator::buildTierInfo(
        const MeshVolume& mesh,
        const PMGroup& pmGroup) {
    assignTiersBasedOnMaxTimeStep(mesh, pmGroup);
    // Grows smallest tier regions for smoothing.
    if (nTiers_ > 1 && growSmallerTiers > 0) {
        growSmallestTierRegions(growSmallerTiers, mesh);
    }
    if (nTiers_ > 1) {
        assignStages(mesh);
    }
    // Ensures that all elem within same tier/stage are consecutive.
    if (nTiers_ > 1) {
        timeTierList_.sortRows_omp(1,2);
        nTiers_ = timeTierList_.maxValInCol(1) + 1;
    }
}

Interval** Integrator::buildTierRange(
        Interval **range,
        const DynMatrix<UInt>& list) {
    UInt nK = list.nRows();
    static const UInt vS = 2;
    const UInt nStages = getNStages();
    // Allocates memory for tierRange.
    if (range != NULL) {
        for (UInt i = 0; i < nTiers_; i++) {
            delete [] range[i];
        }
        delete [] range;
    }
    range = new pair<UInt,UInt>*[nTiers_];
    for (UInt i = 0; i < nTiers_; i++) {
        range[i] = new pair<UInt,UInt>[nStages];
    }
    // Assigns ranges for tier 0.
    if (nTiers_ == 1) {
        for (UInt i = 0; i < nStages; i++) {
            range[0][i].first = 0;
            range[0][i].second = nK;
        }
    } else {
        static const UInt nextKey[2] = {1, 0};
        UInt e2 = list.findFirstOcurrenceInColumns(nextKey,1,vS);
        for (UInt i = 0; i < nStages; i++) {
            range[0][i].first = 0;
            range[0][i].second = e2;
        }
    }
    //
    UInt key[2], nextKey[2];
    for (UInt tier = 1; tier < nTiers_; tier++) {
        for (UInt stage = 0; stage < nStages; stage++) {
            key[0] = tier;
            key[1] = stage;
            UInt e1 = list.findFirstOcurrenceInColumns(key,1,vS);
            UInt e2;
            if (tier+1 != nTiers_ || stage+1 != nStages) {
                nextKey[0] = tier + (stage+1) / nStages;
                nextKey[1] = (stage+1) % nStages;
                e2 = list.findFirstOcurrenceInColumns(nextKey,1,vS);
            } else {
                e2 = nK;
            }
            range[tier][stage].first = e1;
            range[tier][stage].second = e2;
        }
    }
    return range;
}

void Integrator::checkMaterialStabilityForDT(
        const PhysicalModel* mat,
        const Real dt) const {
    const PMVolumeDispersive* disp = dynamic_cast<const PMVolumeDispersive*>(mat);
    if (disp != NULL) {
        for (UInt p = 0; p < disp->getPoleNumber(); p++) {
            Real polePeriod = 1.0 / std::abs(disp->getPole(p));
            if (polePeriod < dt) {
                cerr << endl << "ERROR@Integrator: "
                        << " Contains pole #" << p + 1
                        << " with value " << disp->getPole(p)
                        << " will cause an unstability for dt"
                        << dt << endl;
                mat->printInfo();
            }
        }
    }
}

void Integrator::growSmallestTierRegions(
        const UInt toGrow,
        const MeshVolume& mesh) {
    timeTierList_.sortRows_omp(0,0);
    const UInt nK = mesh.elems().getOf<VolR>().size();
    for (UInt tier = 0; tier < nTiers_-1; tier++) {
        // Creates a list with all the elements belonging to this tier.
        Group<const VolR> elemsInTier;
        elemsInTier.reserve(nK);
        for (UInt k = 0; k < nK; k++) {
            if (timeTierList_(k,1) == tier) {
                ElementId id(timeTierList_(k,0));
                elemsInTier.add(mesh.elems().getId(id));
            }
        }
        Group<const VolR> neigh;
        Group<const VolR> grownElem(elemsInTier);
        for (UInt stage = 0; stage < toGrow; stage++) {
            Group<const VolR> newNeigh = mesh.getAdjacentRegion(grownElem);
            neigh.add(newNeigh);
            grownElem.add(newNeigh);
        }
        for (UInt k = 0; k < neigh.size(); k++) {
            UInt id = neigh(k)->getId().toUInt();
            const UInt row = timeTierList_.findFirstOcurrenceInColumns(&id,0,1);
            timeTierList_(row,1) = tier;
        }
    }
    // Updates number of tiers.
    timeTierList_.sortRows_omp(1,1);
    nTiers_ = timeTierList_(nK-1,1) + 1;
}

void Integrator::assignTiersBasedOnMaxTimeStep(
        const MeshVolume& mesh,
        const PMGroup& pmGroup) {
    // Takes only elements with a defined material.
    GroupElements<const VolR> vol = mesh.elems().getOf<VolR>();
    vol.removeMatId(MatId(0));
    const UInt nK =  vol.size();
    // Computes maximum global timestep (local minimum).
    DynMatrix<Real> dtList(vol.size(), 4);
    mindt = 0.0;
    for (UInt k = 0; k < vol.size(); k++) {
        const PhysicalModel* mat = pmGroup.getId(vol(k)->getMatId());
        Real dt = getMaxTimeStep(vol(k), mat);
        if (mindt > dt || mindt == 0.0) {
            mindt = dt;
        }
        dtList(k,0) = vol(k)->getId().toUInt();
        dtList(k,1) = dt;
        dtList(k,2) = noTier;
        dtList(k,3) = getNStages() - 1;
    }
    if (maxNumOfTiers != 1) {
        dtList.sortRows_omp(1,1);
        Real maxdtList = dtList(nK-1,1);
        Real ratio = getMaxTimeRatio();
        nTiers_ = floor(log(mindt/maxdtList)/log(ratio)) + 1;
        if (maxNumOfTiers > 0 && nTiers_ > maxNumOfTiers) {
            nTiers_ = maxNumOfTiers;
        }
        for (UInt tier = 0; tier < nTiers_; tier++) {
            Real inf = mindt / pow(ratio, Real(tier));
            Real sup;
            if (tier+1 == maxNumOfTiers) {
                sup = numeric_limits<Real>::max();
            } else {
                sup = mindt / pow(ratio, Real(tier + 1));
            }
            for (UInt k = 0; k < nK; k++) {
                if (dtList(k,1) >= inf && dtList(k,1) < sup) {
                    dtList(k,2) = tier;
                }
            }
        }
    } else {
        nTiers_ = 1;
    }
    timeTierList_.copy(dtList.eliminateColumns(1,1));
}

void Integrator::assignStages(const MeshVolume& mesh) {
    if (nTiers_ == 1) {
        return;
    }
    // ----------- Reassigns tiers ------------------------------------
    GroupElements<const VolR> vol = mesh.elems().getOf<VolR>();
    const UInt nK = vol.size();
    timeTierList_.sortRows_omp(0,0);
    const UInt nStages = getNStages();
    for (UInt tier = 0; tier < nTiers_-2; tier++) {
        // Creates a list with all the elements belonging to this tier.
        GroupElements<const VolR> elem;
        elem.reserve(nK);
        bool isInRegion;
        for (UInt k = 0; k < nK; k++) {
            isInRegion = (timeTierList_(k,1) == tier);
            if (isInRegion) {
                elem.add(vol.getId(ElementId(timeTierList_(k,0))));
            }
        }
        GroupElements<const VolR> grownElem = elem;
        GroupElements<const VolR> neigh;
        for (UInt stage = 0; stage < (nStages * growStages); stage++) {
            GroupElements<const VolR> newNeigh = mesh.getAdjacentRegion(grownElem);
            neigh.add(newNeigh);
            grownElem.add(newNeigh);
        }
        for (UInt k = 0; k < neigh.size(); k++) {
            UInt id = neigh(k)->getId().toUInt();
            UInt row = timeTierList_.findFirstOcurrenceInColumns(&id,0,1);
            if (timeTierList_(row,1) > tier + 1) {
                timeTierList_(row,1) = tier + 1;
            }
        }
    }
    // ----------- Assigns stages -------------------------------------
    for (UInt tier = 0; tier < nTiers_-1; tier++) {
        GroupElements<const VolR> elem;
        elem.reserve(nK);
        bool isInRegion;
        for (UInt k = 0; k < nK; k++) {
            isInRegion = (timeTierList_(k,1) == tier);
            if (isInRegion) {
                elem.add(vol.getId(ElementId(timeTierList_(k,0))));
            }
        }
        GroupElements<const VolR> grownElem = elem;
        GroupElements<const VolR> neigh;
        for (UInt stage = 0; stage < nStages; stage++) {
            for (UInt times = 0; times < growStages; times++) {
                GroupElements<const VolR> newNeigh = mesh.getAdjacentRegion(grownElem);
                grownElem.add(newNeigh);
                for (UInt k = 0; k < newNeigh.size(); k++) {
                    UInt id = newNeigh(k)->getId().toUInt();
                    UInt row =
                            timeTierList_.findFirstOcurrenceInColumns(&id,0,1);
                    if (row < nK) {
                        if (timeTierList_(row,1) > tier) {
                            timeTierList_(row,2) = stage;
                        }
                    }
                }
            }
        }
    }
}

vector<pair<UInt, UInt>> Integrator::getIdPartitionVector(
        const vector<vector<ElementId>>& pId) const {
    vector<pair<UInt,UInt>> res;
    UInt nPart = pId.size();
    UInt gSize = 0;
    for (UInt p = 0; p < nPart; p++) {
        gSize += pId[p].size();
    }
    res.reserve(gSize);
    pair<UInt,UInt> idPart;
    for (UInt p = 0; p < nPart; p++) {
        for (UInt k = 0; k < pId[p].size(); k++) {
            idPart.first = pId[p][k].toUInt();
            idPart.second = p;
            res.push_back(idPart);
        }
    }
    return res;
}

