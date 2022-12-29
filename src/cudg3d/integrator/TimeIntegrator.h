#pragma once

#include "cudg3d/dg/Evolution.h"

namespace SEMBA::dgtd::integrator {

typedef std::pair<size_t,size_t> Interval;

class TimeIntegrator {
public:
    struct Options {
        enum class Type {
            lserk4, lf2
        };
        Type timeIntegrator{ Type::lserk4 };

        Math::Real finalTime = 30e-9;

        bool useLTS{ false };
        std::size_t growSmallerTiers = 0;
        std::size_t maxNumberOfTiers = 0;
        bool useMaxStageSizeForLTS = false;
    };
    
    TimeIntegrator(const Options&);
    virtual ~TimeIntegrator() = default;
    
    //virtual void timeIntegrate(const Math::Real timer) const = 0;
    //void setSolver(DG* solver);
    //Math::Real getMaxDT() const;
    //Math::Real getMinDT() const;
    //size_t getNTiers() const;
    //size_t getNPartitions() const;
    //vector<vector<Geometry::ElemId>> getTiersIds() const;
    //vector<vector<Geometry::ElemId>> getStagesIds() const;
    //vector<vector<Geometry::ElemId>> getPartitionsIds() const;
    //Interval getRange(const size_t tier, const size_t stage) const;
    //vector<pair<Geometry::ElemId,Math::Int>> getComputationalWeights(
    //        const Mesh::Volume* msh) const;
    //void partitionate(const Mesh::Volume* mesh, Communications::Comm* comm);

private:
    Options opts_;
//    static const size_t noTier = 0;
//    static const size_t noStage = 0;
//    static const size_t growStages = 1;
//    size_t nTiers_;
//    pair<size_t,size_t> **tierRange_;
//    vector<vector<Geometry::ElemId>> partIds_;
//    DG* solver;
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
};

}
