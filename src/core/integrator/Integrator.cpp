//// OpenSEMBA
//// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
////                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
////                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
////                    Daniel Mateos Romero            (damarro@semba.guru)
////
//// This file is part of OpenSEMBA.
////
//// OpenSEMBA is free software: you can redistribute it and/or modify it under
//// the terms of the GNU Lesser General Public License as published by the Free
//// Software Foundation, either version 3 of the License, or (at your option)
//// any later version.
////
//// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
//// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//// details.
////
//// You should have received a copy of the GNU Lesser General Public License
//// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
///*
// * Integrator.cpp
// *
// *  Created on: Feb 21, 2013
// *      Author: luis
// */
//
//#include "Integrator.h"
//
//namespace SEMBA {
//namespace Cudg3d {
//namespace Integrator {
//
//Integrator::Integrator() {
//    mindt = 0.0;
//    nTiers_ = 0;
//    tierRange_ = NULL;
//    doLTS = true;
//    growSmallerTiers = 0;
//    maxNumOfTiers = 0;
//    cfl_ = 0.0;
//    solver = NULL;
//}
//
//Integrator::~Integrator() {
//    for (size_t i = 0; i < nTiers_; i++) {
//        delete tierRange_[i];
//    }
//}
//
//void Integrator::setSolver(DG* solver_) {
//    solver = solver_;
//}
//
//size_t Integrator::getNTiers() const {
//    return nTiers_;
//}
//
//Math::Real Integrator::getMaxDT() const {
//    if (doLTS) {
//        return (mindt / pow(getMaxTimeRatio(), Math::Real(nTiers_-1)));
//    } else {
//        return mindt;
//    }
//}
//
//Math::Real Integrator::getMinDT() const {
//    return mindt;
//}
//
//vector<vector<Geometry::ElemId>> Integrator::getTiersIds() const {
//    vector<vector<Geometry::ElemId>> res;
//    for (size_t tier = 0; tier < nTiers_; tier++) {
//        vector<Geometry::ElemId> tierIds = getIdsOfTier(tier);
//        res.push_back(tierIds);
//    }
//    return res;
//}
//
//vector<vector<Geometry::ElemId>> Integrator::getStagesIds() const {
//    vector<vector<Geometry::ElemId> > res;
//    for (size_t stage = 0; stage < getNStages(); stage++) {
//        vector<Geometry::ElemId> stageIds = getIdsOfStage(stage);
//        res.push_back(stageIds);
//    }
//    return res;
//}
//
//vector<vector<Geometry::ElemId>> Integrator::getPartitionsIds() const {
//    //assert(partIds.size() != 0);
//    return partIds_;
//}
//
//vector<pair<Geometry::ElemId,Math::Int>> Integrator::getComputationalWeights(
//        const Mesh::Volume* msh) const {
//    const Math::Int curlFlops = 1;
//    const Math::Int fluxFlops = 0;
//    Math::Int flops = curlFlops + fluxFlops;
//    Geometry::VolRGroup physVol = msh->elems();
//    physVol.removeMatId(MatId(0));
//    const size_t nK = physVol.sizeOf<Geometry::VolR>();
//    vector<pair<Geometry::ElemId,Math::Int>> idWgt;
//    idWgt.reserve(nK);
//    for (size_t e = 0; e < nK; e++) {
//        pair<Geometry::ElemId,Math::Int> aux;
//        Geometry::ElemId id = getIdOfGlobalRelPos(e);
//        aux.first = id;
//        aux.second = getNumOfIterationsPerBigTimeStep(e) * flops;
//        idWgt.push_back(aux);
//    }
//    return idWgt;
//}
//
//Interval Integrator::getRange(const size_t tier, const size_t stage) const {
//    assert(tier < nTiers_);
//    assert(stage < getNStages());
//    assert(tierRange_ != NULL);
//    return tierRange_[tier][stage];
//}
//
//size_t Integrator::getNPartitions() const {
//    return partIds_.size();
//}
//
//void Integrator::partitionate(
//        const Mesh::Volume* msh,
//        Communications::Comm* comm) {
//    cout << " - Getting computational weights... " << flush;
//    vector<pair<Geometry::ElemId,Math::Int> > idWgt =
//            getComputationalWeights(msh);
//    cout << "OK" << endl;
//    cout << " - Obtaining partition ids... " << flush;
//    vector<vector <Geometry::ElemId>> partId =
//            msh->getPartitionsIds(comm->getNumberOfTasks(), idWgt);
//    cout << "OK" << endl;
//    cout << " - Setting partition sizes... " << flush;
//    comm->setPartitionSizes(partId);
//    cout << "OK" << endl;
//    cout << " - Reordering partitions... " << flush;
//    reorder(partId, comm->getLocalOffset(), comm->getLocalSize());
//    cout << "OK" << endl;
//}
//
//void Integrator::printInfo() const {
//    cout << "--- SolverInfo ---" << endl;
//    cout << "Min. time step: " << mindt*1E12 << " [ps]" << endl;
//    cout << "Max. time step: " << getMaxDT()*1E12 << " [ps]" << endl;
//    cout << "Number of tiers: " << nTiers_ << endl;
//    if (nTiers_ > 1) {
//        for (size_t i = 0; i < nTiers_; i++) {
//            cout << "# of Cells in tier " << i << ": "
//                    << getNumberOfCellsInTier(i) << endl;
//            for (size_t j = 0; j < getNStages(); j++) {
//                cout << "--> Range stage #" << j << ": "
//                        << getRange(i,j).first << ", "
//                        << getRange(i,j).second << endl;
//            }
//        }
//    }
//}
//
//void Integrator::init(
//        const Mesh::Volume& mesh,
//        const PMGroup& pmGroup,
//        const OptionsSolverDGTD* arg) {
//    growSmallerTiers = arg->getGrowSmallerTiers();
//    maxNumOfTiers = arg->getMaxNumberOfTiers();
//    doLTS = arg->isUseLTS();
//    if (!doLTS) {
//        maxNumOfTiers = 1;
//    }
//    cout << "- Building tier info... " << flush;
//    buildTierInfo(mesh, pmGroup);
//    cout << " OK" << endl;
//    cout << "- Building Relative Positions of Ids... " << flush;
//    buildRelPosOfIds(timeTierList_);
//    cout << " OK" << endl;
//    cout << "- Building Tier Range... " << flush;
//    tierRange_ = buildTierRange(tierRange_, timeTierList_);
//    cout <<  "OK" << endl;
//}
//
//Math::Real Integrator::getMaxTimeStep(
//        const Geometry::VolR* vol,
//        const PhysicalModel* mat) const {
//    // Returns the maximum time step allowed for this cell.
//    Math::Real fS1 = 0.0;
//    for (size_t f = 0; f < vol->numberOfFaces(); f++) {
//        Math::Real area = vol->getAreaOfFace(f);
//        Math::Real volume = vol->getVolume();
//        Math::Real fS2 = area / volume;
//        if (fS2 > fS1) {
//            fS1 = fS2;
//        }
//    }
//    Math::Real dt = (1.0 / Math::Constants::c0) *
//            (1.0 / Math::Real(fS1 * ORDER_N * ORDER_N));
//    // Checks case of electrical dispersive materials.
//    checkMaterialStabilityForDT(mat, dt);
//    return (dt * cfl_);
//}
//
//size_t Integrator::getNumberOfCellsInTier(const size_t tier) const {
//    assert(tier < nTiers_);
//    assert(tierRange_ != NULL);
//    size_t first = getRange(tier, 0).first;
//    size_t last = getRange(tier, getNStages() - 1).second;
//    size_t res = last - first;
//    return res;
//}
//
//vector<Geometry::ElemId> Integrator::getIdsOfTier(size_t tier) const {
//    assert(tier < nTiers_);
//    const size_t nK = timeTierList_.nRows();
//    vector<Geometry::ElemId> res;
//    res.reserve(nK);
//    for (size_t i = 0; i < nK; i++) {
//        if (timeTierList_(i,1) == tier) {
//            res.push_back(Geometry::ElemId(timeTierList_(i,0)));
//        }
//    }
//    return res;
//}
//
//vector<Geometry::ElemId> Integrator::getIdsOfStage(size_t stage) const {
//    assert(stage < getNStages());
//    size_t nK = timeTierList_.nRows();
//    vector<Geometry::ElemId> res;
//    res.reserve(nK);
//    for (size_t i = 0; i < nK; i++) {
//        if (timeTierList_(i,2) == stage) {
//            res.push_back(Geometry::ElemId(timeTierList_(i,0)));
//        }
//    }
//    return res;
//}
//
//void Integrator::reorder(
//        const vector<vector<Geometry::ElemId> >& partitionId,
//        const size_t localOffset,
//        const size_t localSize) {
//    partIds_ = partitionId;
//    // Reorders timeTierList according to partitions.
//    reorderTimeTierList(partitionId);
//    // Now aux stores the final ordering considering partitions.
//    buildRelPosOfIds(timeTierList_);
//    // Builds time tier range
//    Math::Matrix::Dynamic   <size_t> reducedList(localSize,3);
//    for (size_t i = 0; i < localSize; i++) {
//        for (size_t j = 0; j < 3; j++) {
//            reducedList(i,j) = timeTierList_(i + localOffset, j);
//        }
//    }
//    tierRange_ = buildTierRange(tierRange_, reducedList);
//}
//
//void Integrator::reorderTimeTierList(
//        const vector<vector<Geometry::ElemId>>& partitionId) {
//    size_t nK = timeTierList_.nRows();
//    Math::Matrix::Dynamic<size_t> aux(nK, 5); // relPos - Ids - Part - Tier - Stage
//    for (size_t k = 0; k < nK; k++) {
//        aux(k, 0) = k;
//        aux(k, 1) = timeTierList_(k, 0);
//        aux(k, 2) = 0; // Temporary assignation of partition.
//        aux(k, 3) = timeTierList_(k, 1); // Tier assignation.
//        aux(k, 4) = timeTierList_(k, 2);
//    }
//    aux.sortRows_omp(1, 1);
//    size_t initId = aux(0, 1);
//    size_t nParts = partitionId.size();
//    for (size_t p = 0; p < nParts; p++) {
//        for (size_t i = 0; i < partitionId[p].size(); i++) {
//            size_t id = partitionId[p][i];
//            aux(id - initId, 2) = p;
//        }
//    }
//    aux.sortRows_omp(2, 4);
//    for (size_t i = 0; i < nK; i++) {
//        timeTierList_(i, 0) = aux(i, 1);
//        timeTierList_(i, 1) = aux(i, 3);
//        timeTierList_(i, 2) = aux(i, 4);
//    }
//}
//
//void Integrator::buildTierInfo(
//        const Mesh::Volume& mesh,
//        const PMGroup& pmGroup) {
//    assignTiersBasedOnMaxTimeStep(mesh, pmGroup);
//    // Grows smallest tier regions for smoothing.
//    if (nTiers_ > 1 && growSmallerTiers > 0) {
//        growSmallestTierRegions(growSmallerTiers, mesh);
//    }
//    if (nTiers_ > 1) {
//        assignStages(mesh);
//    }
//    // Ensures that all elem within same tier/stage are consecutive.
//    if (nTiers_ > 1) {
//        timeTierList_.sortRows_omp(1,2);
//        nTiers_ = timeTierList_.maxValInCol(1) + 1;
//    }
//}
//
//Interval** Integrator::buildTierRange(
//        Interval **range,
//        const Math::Matrix::Dynamic<size_t>& list) {
//    size_t nK = list.nRows();
//    static const size_t vS = 2;
//    const size_t nStages = getNStages();
//    // Allocates memory for tierRange.
//    if (range != NULL) {
//        for (size_t i = 0; i < nTiers_; i++) {
//            delete [] range[i];
//        }
//        delete [] range;
//    }
//    range = new pair<size_t,size_t>*[nTiers_];
//    for (size_t i = 0; i < nTiers_; i++) {
//        range[i] = new pair<size_t,size_t>[nStages];
//    }
//    // Assigns ranges for tier 0.
//    if (nTiers_ == 1) {
//        for (size_t i = 0; i < nStages; i++) {
//            range[0][i].first = 0;
//            range[0][i].second = nK;
//        }
//    } else {
//        static const size_t nextKey[2] = {1, 0};
//        size_t e2 = list.findFirstOcurrenceInColumns(nextKey,1,vS);
//        for (size_t i = 0; i < nStages; i++) {
//            range[0][i].first = 0;
//            range[0][i].second = e2;
//        }
//    }
//    //
//    size_t key[2], nextKey[2];
//    for (size_t tier = 1; tier < nTiers_; tier++) {
//        for (size_t stage = 0; stage < nStages; stage++) {
//            key[0] = tier;
//            key[1] = stage;
//            size_t e1 = list.findFirstOcurrenceInColumns(key,1,vS);
//            size_t e2;
//            if (tier+1 != nTiers_ || stage+1 != nStages) {
//                nextKey[0] = tier + (stage+1) / nStages;
//                nextKey[1] = (stage+1) % nStages;
//                e2 = list.findFirstOcurrenceInColumns(nextKey,1,vS);
//            } else {
//                e2 = nK;
//            }
//            range[tier][stage].first = e1;
//            range[tier][stage].second = e2;
//        }
//    }
//    return range;
//}
//
//void Integrator::checkMaterialStabilityForDT(
//        const PhysicalModel::PhysicalModel* mat,
//        const Math::Real dt) const {
//    if (!mat->is<PhysicalModel::Volume::Dispersive>()) {
//        return;
//    }
//    const PhysicalModel::Volume::Dispersive* disp =
//            mat->castTo<PhysicalModel::Volume::Dispersive>();
//    for (size_t p = 0; p < disp->getPoleNumber(); p++) {
//        Math::Real polePeriod = 1.0 / std::abs(disp->getPole(p));
//        if (polePeriod < dt) {
//            cerr << endl << "ERROR@Integrator: "
//                    << " Contains pole #" << p + 1
//                    << " with value " << disp->getPole(p)
//                    << " will cause an unstability for dt"
//                    << dt << endl;
//            mat->printInfo();
//            throw logic_error("Invalid DT for pole.");
//        }
//    }
//}
//
//void Integrator::growSmallestTierRegions(
//        const size_t toGrow,
//        const Mesh::Volume& mesh) {
//    timeTierList_.sortRows_omp(0,0);
//    const size_t nK = mesh.elems().sizeOf<Geometry::VolR>();
//    for (size_t tier = 0; tier < nTiers_-1; tier++) {
//        // Creates a list with all the elements belonging to this tier.
//        Group<const Geometry::VolR> elemsInTier;
//        elemsInTier.reserve(nK);
//        for (size_t k = 0; k < nK; k++) {
//            if (timeTierList_(k,1) == tier) {
//                Geometry::ElemId id(timeTierList_(k,0));
//                elemsInTier.add(mesh.elems().getId(id));
//            }
//        }
//        Group<const Geometry::VolR> neigh;
//        Group<const Geometry::VolR> grownElem(elemsInTier);
//        for (size_t stage = 0; stage < toGrow; stage++) {
//            Group<const Geometry::VolR> newNeigh =
//                    mesh.getAdjacentRegion(grownElem);
//            neigh.add(newNeigh);
//            grownElem.add(newNeigh);
//        }
//        for (size_t k = 0; k < neigh.size(); k++) {
//            size_t id = neigh(k)->getId().tosize_t();
//            const size_t row = timeTierList_.findFirstOcurrenceInColumns(&id,0,1);
//            timeTierList_(row,1) = tier;
//        }
//    }
//    // Updates number of tiers.
//    timeTierList_.sortRows_omp(1,1);
//    nTiers_ = timeTierList_(nK-1,1) + 1;
//}
//
//void Integrator::assignTiersBasedOnMaxTimeStep(
//        const Mesh::Volume& mesh,
//        const PMGroup& pmGroup) {
//    // Takes only elements with a defined material.
//    Geometry::VolRGroup vol = mesh.elems().getOf<VolR>();
//    vol.removeMatId(MatId(0));
//    const size_t nK =  vol.size();
//    // Computes maximum global timestep (local minimum).
//    Math::Matrix::Dynamic<Math::Real> dtList(vol.size(), 4);
//    mindt = 0.0;
//    for (size_t k = 0; k < vol.size(); k++) {
//        const PhysicalModel* mat = pmGroup.getId(vol(k)->getMatId());
//        Math::Real dt = getMaxTimeStep(vol(k), mat);
//        if (mindt > dt || mindt == 0.0) {
//            mindt = dt;
//        }
//        dtList(k,0) = vol(k)->getId().tosize_t();
//        dtList(k,1) = dt;
//        dtList(k,2) = noTier;
//        dtList(k,3) = getNStages() - 1;
//    }
//    if (maxNumOfTiers != 1) {
//        dtList.sortRows_omp(1,1);
//        Math::Real maxdtList = dtList(nK-1,1);
//        Math::Real ratio = getMaxTimeRatio();
//        nTiers_ = floor(log(mindt/maxdtList)/log(ratio)) + 1;
//        if (maxNumOfTiers > 0 && nTiers_ > maxNumOfTiers) {
//            nTiers_ = maxNumOfTiers;
//        }
//        for (size_t tier = 0; tier < nTiers_; tier++) {
//            Math::Real inf = mindt / pow(ratio, Math::Real(tier));
//            Math::Real sup;
//            if (tier+1 == maxNumOfTiers) {
//                sup = numeric_limits<Math::Real>::max();
//            } else {
//                sup = mindt / pow(ratio, Math::Real(tier + 1));
//            }
//            for (size_t k = 0; k < nK; k++) {
//                if (dtList(k,1) >= inf && dtList(k,1) < sup) {
//                    dtList(k,2) = tier;
//                }
//            }
//        }
//    } else {
//        nTiers_ = 1;
//    }
//    timeTierList_.copy(dtList.eliminateColumns(1,1));
//}
//
//void Integrator::assignStages(const Mesh::Volume& mesh) {
//    if (nTiers_ == 1) {
//        return;
//    }
//    // ----------- Reassigns tiers ------------------------------------
//    Geometry::VolRGroup vol = mesh.elems().getOf<Geometry::VolR>();
//    const size_t nK = vol.size();
//    timeTierList_.sortRows_omp(0,0);
//    const size_t nStages = getNStages();
//    for (size_t tier = 0; tier < nTiers_-2; tier++) {
//        // Creates a list with all the elements belonging to this tier.
//        Geometry::VolRGroup elem;
//        elem.reserve(nK);
//        bool isInRegion;
//        for (size_t k = 0; k < nK; k++) {
//            isInRegion = (timeTierList_(k,1) == tier);
//            if (isInRegion) {
//                elem.add(vol.getId(Geometry::ElemId(timeTierList_(k,0))));
//            }
//        }
//        Geometry::VolRGroup grownElem = elem;
//        Geometry::VolRGroup neigh;
//        for (size_t stage = 0; stage < (nStages * growStages); stage++) {
//            Geometry::VolRGroup newNeigh = mesh.getAdjacentRegion(grownElem);
//            neigh.add(newNeigh);
//            grownElem.add(newNeigh);
//        }
//        for (size_t k = 0; k < neigh.size(); k++) {
//            size_t id = neigh(k)->getId();
//            size_t row = timeTierList_.findFirstOcurrenceInColumns(&id,0,1);
//            if (timeTierList_(row,1) > tier + 1) {
//                timeTierList_(row,1) = tier + 1;
//            }
//        }
//    }
//    // ----------- Assigns stages -------------------------------------
//    for (size_t tier = 0; tier < nTiers_-1; tier++) {
//        Geometry::VolRGroup elem;
//        elem.reserve(nK);
//        bool isInRegion;
//        for (size_t k = 0; k < nK; k++) {
//            isInRegion = (timeTierList_(k,1) == tier);
//            if (isInRegion) {
//                elem.add(vol.getId(Geometry::ElemId(timeTierList_(k,0))));
//            }
//        }
//        Geometry::VolRGroup grownElem = elem;
//        Geometry::VolRGroup neigh;
//        for (size_t stage = 0; stage < nStages; stage++) {
//            for (size_t times = 0; times < growStages; times++) {
//                Geometry::VolRGroup newNeigh = mesh.getAdjacentRegion(grownElem);
//                grownElem.add(newNeigh);
//                for (size_t k = 0; k < newNeigh.size(); k++) {
//                    size_t id = newNeigh(k)->getId().tosize_t();
//                    size_t row =
//                            timeTierList_.findFirstOcurrenceInColumns(&id,0,1);
//                    if (row < nK) {
//                        if (timeTierList_(row,1) > tier) {
//                            timeTierList_(row,2) = stage;
//                        }
//                    }
//                }
//            }
//        }
//    }
//}
//
//vector<pair<size_t, size_t>> Integrator::getIdPartitionVector(
//        const vector<vector<Geometry::ElemId>>& pId) const {
//    vector<pair<size_t,size_t>> res;
//    size_t nPart = pId.size();
//    size_t gSize = 0;
//    for (size_t p = 0; p < nPart; p++) {
//        gSize += pId[p].size();
//    }
//    res.reserve(gSize);
//    pair<size_t,size_t> idPart;
//    for (size_t p = 0; p < nPart; p++) {
//        for (size_t k = 0; k < pId[p].size(); k++) {
//            idPart.first = pId[p][k].tosize_t();
//            idPart.second = p;
//            res.push_back(idPart);
//        }
//    }
//    return res;
//}
//
//}
//}
//}
