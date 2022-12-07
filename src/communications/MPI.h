#ifndef COMMUNICATIONSMPI_H_
#define COMMUNICATIONSMPI_H_

#include "MPI.h"
#include "Communications.h"


namespace SEMBA {
namespace Cudg3d {

struct Field_s {
	Math::Int rp;
	Math::Real Ex, Ey, Ez, Hx, Hy, Hz;
};

//class MPI : public Communications {
//public:
//	MPI();
//	virtual ~MPI();
//	Math::Int
//	 getNumberOfTasks() const;
//	Math::Int
//	 getTask() const;
//	size_t
//	 getLocalOffset() const;
//	size_t
//	 getLocalSize() const;
//	Math::Int
//	 getNumOfTasksOnThisHost() const;
//	void
//	 setPartitionSizes(
//	  const vector<vector<size_t> >& partId);
//	void
//	 gatherFieldsMaster(
// 	  Math::Real* Ex, Math::Real* Ey, Math::Real* Ez,
//	  Math::Real* Hx, Math::Real* Hy, Math::Real* Hz,
//	  const Math::Real* lEx, const Math::Real* lEy, const Math::Real* lEz,
//	  const Math::Real* lHx, const Math::Real* lHy, const Math::Real* lHz) const;
//	void
//	 gatherFieldsSlave(
// 	  const Math::Real* Ex, const Math::Real* Ey, const Math::Real* Ez,
//	  const Math::Real* Hx, const Math::Real* Hy, const Math::Real* Hz) const;
//	bool
//	 isMaster() const;
//	void
//	 syncNeighbourFields(
//	  Math::Real* nEx, Math::Real* nEy, Math::Real* nEz,
//	  Math::Real* nHx, Math::Real* nHy, Math::Real* nHz,
//	  const Math::Real* Ex, const Math::Real* Ey, const Math::Real* Ez,
//	  const Math::Real* Hx, const Math::Real* Hy, const Math::Real* Hz) const;
//	void
//	 initNeighbourFields(const vector<size_t>& nIds);
//	Math::Real
//	 reduceToGlobalMinimum(Math::Real val) const;
//	void
//	 prMath::IntInfo() const;
//private:
//	static const Math::Int master = 0;
//	static const size_t np = ((ORDER_N+1)*(ORDER_N+2)*(ORDER_N+3)/6);
//	static const bool weightPartitions = true;
//	Math::Int nTasks;
//	Math::Int nTasksInHost;
//	Math::Int task;
//	Math::Int len;
//	Math::Int* pSize;
//	Math::Int* pOffset;
//	Math::Int* fSize;
//	Math::Int* fOffset;
//	size_t* neighElemRP;
//	MPI_Comm world;
//	MPI_Datatype MPIField;
//	size_t nPartitions;
//	Math::Int* nNeighId;
//	size_t** neighId;
//	Math::Int **neighFSize;
//	Math::Int **neighFOffset;
//	Math::Int **nDofSize;
//	Math::Int **nDofOffset;
//	Math::Int*
//	 getDofSizes() const;
//	Math::Int*
//	 getDofOffsets() const;
//	size_t
//	 getLocalDofSize() const;
//	size_t
//	 getGlobalDofSize() const;
//	void
//	 packFields(
//	  Field_s *field,
//	  const Math::Real *Ex, const Math::Real *Ey, const Math::Real *Ez,
//	  const Math::Real *Hx, const Math::Real *Hy, const Math::Real *Hz,
//	  const size_t fSize) const;
//	void
//	 unpackFields(
//	  Math::Real *Ex, Math::Real *Ey, Math::Real *Ez,
//	  Math::Real *Hx, Math::Real *Hy, Math::Real *Hz,
//	  const Field_s *field, const size_t fSize) const;
//	void
//	 commitFieldStruct();
//	Math::Int*
//	 getFieldSizes() const;
//	Math::Int*
//	 getFieldOffsets() const;
//	size_t
//	 getLocalFieldSize() const;
//	Math::Int
//	 getTaskOfId(const size_t id) const;
//	bool
//	 checkNNeighCoherence(Math::Int* nneigh) const;
//	bool
//	 checkVectorsAreEqual(
//	  const size_t vSize,
//	  const size_t* v,
//	  const vector<size_t>& nIds) const;
//	bool
//	 checkNeighFSizes() const;
//	Math::Int
//	 countTasksInLocalHost() const;
//	vector<size_t>
//	 getThreadsOfTasks() const;
//};

}
}

#endif /* COMMUNICATIONSMPI_H_ */
