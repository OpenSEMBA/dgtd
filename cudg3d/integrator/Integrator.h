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
// * SolverInfo.h
// *
// *  Created on: Feb 21, 2013
// *      Author: luis
// */
//
//#ifndef INTEGRATOR_H_
//#define INTEGRATOR_H_
//
//#include <vector>
//#include <utility>
//#include <limits>
//
//using namespace std;
//
//#include "core/Comm.h"
//#include "core/mesh/Volume.h"
//#include "physicalModel/PhysicalModel.h"
//#include "options/OptionsSolverDGTD.h"
//#include "DG.h"
//
//namespace SEMBA {
//namespace Cudg3d {
//namespace Integrator {
//
//typedef pair<size_t,size_t> Interval;
//
//class Integrator {
//public:
//    Integrator();
//    virtual ~Integrator();
//    virtual void timeIntegrate(const Math::Real timer) const = 0;
//    void setSolver(DG* solver);
//    Math::Real getMaxDT() const;
//    Math::Real getMinDT() const;
//    size_t getNTiers() const;
//    size_t getNPartitions() const;
//    vector<vector<Geometry::ElemId>> getTiersIds() const;
//    vector<vector<Geometry::ElemId>> getStagesIds() const;
//    vector<vector<Geometry::ElemId>> getPartitionsIds() const;
//    Interval getRange(const size_t tier, const size_t stage) const;
//    vector<pair<Geometry::ElemId,Math::Int>> getComputationalWeights(
//            const Mesh::Volume* msh) const;
//    void partitionate(const Mesh::Volume* mesh, Communications::Comm* comm);
//    void printInfo() const;
//protected:
//    Math::Real cfl_;
//    DG* solver;
//    bool doLTS;
//    Math::Matrix::Dynamic<size_t> timeTierList_; // Id - Tier - Stage
//    Math::Real mindt;
//    void init(
//            const Mesh::Volume& mesh,
//            const PMGroup& pmGroup,
//            const Options* arg);
//    size_t getNumberOfCellsInTier(const size_t tier) const;
//    virtual size_t getNumOfIterationsPerBigTimeStep(
//            const size_t e) const = 0;
//    virtual size_t getNStages() const = 0;
//    virtual Math::Real getMaxTimeStep(
//            const Geometry::VolR* tet,
//            const PhysicalModel* mat) const;
//    virtual Math::Real getMaxTimeRatio() const = 0;
//    vector<Geometry::ElemId> getIdsOfTier(const size_t tier) const;
//    vector<Geometry::ElemId> getIdsOfStage(const size_t stage) const;
//private:
//    static const size_t noTier = 0;
//    static const size_t noStage = 0;
//    static const size_t growStages = 1;
//    size_t growSmallerTiers;
//    size_t maxNumOfTiers;
//    size_t nTiers_;
//    pair<size_t,size_t> **tierRange_;
//    vector<vector<Geometry::ElemId>> partIds_;
//    void reorder(
//            const vector<vector<Geometry::ElemId>>& partitionsIds_,
//            const size_t localOffset,
//            const size_t localSize);
//    void buildTierInfo(
//            const Mesh::Volume& mesh,
//            const PMGroup& pmGroup);
//    virtual void checkMaterialStabilityForDT(
//            const PhysicalModel::PhysicalModel* mat,
//            const Math::Real dt) const;
//    void assignTiersBasedOnMaxTimeStep(
//            const Mesh::Volume& mesh,
//            const PMGroup& pmGroup);
//    pair<size_t,size_t>** buildTierRange(
//            pair<size_t,size_t> **range,
//            const Math::Matrix::Dynamic<size_t>& list);
//    void growSmallestTierRegions(
//            const size_t numToGrow,
//            const Mesh::Volume& mesh);
//    vector<pair<size_t, size_t> > getIdPartitionVector(
//            const vector<vector<Geometry::ElemId> >& pId) const;
//    void assignStages(const Mesh::Volume& mesh);
//    void reorderTimeTierList(const vector<vector<Geometry::ElemId>>& partitionId);
//};
//
//
//}
//}
//}
//
//
//#endif /* SOLVERINFO_H_ */
