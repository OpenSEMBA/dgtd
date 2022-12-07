#pragma once

#include <iostream>
#include <assert.h>
#include <vector>

namespace SEMBA {
namespace Cudg3d {
namespace Communications {

class Communications {
public:
    virtual ~Communications() = default;
//    virtual Math::Int getNumberOfTasks() const = 0;
//    virtual Math::Int getTask() const = 0;
//    virtual bool isMaster() const = 0;
//    virtual Math::Int getNumOfTasksOnThisHost() const = 0;
//    virtual size_t getLocalOffset() const = 0;
//    virtual size_t getLocalSize() const = 0;
//    virtual void setPartitionSizes(const vector<vector<Geometry::ElemId>>& partId) = 0;
//    virtual void gatherFieldsMaster(
//            Math::FieldR3& elec, Math::Field<Math::Real,3>& magn,
//            const Math::FieldR3& localElec, const Math::FieldR3& localMagn) const = 0;
//    virtual void gatherFieldsSlave(
//            const Math::FieldR3& electric, const Math::FieldR3& magnetic) const = 0;
//    virtual void syncNeighbourFields(
//            Math::Real* nEx, Math::Real* nEy, Math::Real* nEz,
//            Math::Real* nHx, Math::Real* nHy, Math::Real* nHz,
//            const Math::Real* Ex, const Math::Real* Ey, const Math::Real* Ez,
//            const Math::Real* Hx, const Math::Real* Hy, const Math::Real* Hz) const = 0;
//    virtual void initNeighbourFields(const vector<ElemId>& nIds) = 0;
//    virtual Math::Real reduceToGlobalMinimum(Math::Real val) const = 0;
//    virtual void printInfo() const = 0;
};

}
}
}