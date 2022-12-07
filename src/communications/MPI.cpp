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
//
//
//#include "MPI.h"
//
//namespace SEMBA {
//namespace Cudg3d {
//
//MPI::MPI() {
//	world = MPI_COMM_WORLD;
//	Math::Int rc = MPI_Init(NULL, NULL);
//	if (rc != MPI_SUCCESS) {
//		cout << "Error starting MPI program. Terminating." << endl;
//		MPI_Abort(world, rc);
//    }
//	MPI_Comm_size(world,&nTasks);
//	MPI_Comm_rank(world,&task);
//	commitFieldStruct();
//	pSize = NULL;
//	pOffset = NULL;
//	nPartitions = 1;
//	nNeighId = NULL;
//	neighId = NULL;
//	neighFSize = NULL;
//	neighFOffset = NULL;
//	nDofSize = NULL;
//	nDofOffset = NULL;
//	nTasksInHost = countTasksInLocalHost();
//	len = 0;
//}
//
//MPI::~MPI() {
//	MPI_Finalize();
//}
//
//Math::Int MPI::getNumberOfTasks() const {
//	return nTasks;
//}
//
//Math::Int MPI::getTask() const {
//	return task;
//}
//
//bool MPI::isMaster() const {
//	if (task == master) {
//		return true;
//	}
//	return false;
//}
//
//void MPI::gatherFieldsMaster(
// Math::Real* Ex, Math::Real* Ey, Math::Real* Ez,
// Math::Real* Hx, Math::Real* Hy, Math::Real* Hz,
// const Math::Real* lEx, const Math::Real* lEy, const Math::Real* lEz,
// const Math::Real* lHx, const Math::Real* lHy, const Math::Real* lHz) const {
//	assert(isMaster());
//	// Packs fields in field struct.
//	size_t lFieldSize = getLocalFieldSize();
//	Field_s *lField;
//	lField = new Field_s[lFieldSize];
//	packFields(lField, lEx,lEy,lEz,lHx,lHy,lHz, lFieldSize);
//	// Allocates space for general field.
//	size_t gFieldSize = getGlobalSize() * np;
//	Field_s *gField;
//	gField = new Field_s[gFieldSize];
//	// Gathers fields in master.
//	MPI_Barrier(world);
//	MPI_Gatherv(lField, lFieldSize, MPIField, gField,
//	 getFieldSizes(), getFieldOffsets(), MPIField, master, world);
//	delete lField;
//	// Unpacks in general fields variables.
//	unpackFields(Ex,Ey,Ez,Hx,Hy,Hz,gField, gFieldSize);
//	delete gField;
//}
//
//void
//MPI::gatherFieldsSlave(
// const Math::Real* lEx, const Math::Real* lEy, const Math::Real* lEz,
// const Math::Real* lHx, const Math::Real* lHy, const Math::Real* lHz) const {
//	assert(!isMaster());
//	// Packs fields in field struct.
//	size_t lFSize = getLocalFieldSize();
//	Field_s * lField;
//	lField = new Field_s[lFSize];
//	packFields(lField, lEx,lEy,lEz,lHx,lHy,lHz, lFSize);
//	// Gathers fields in master.
//	MPI_Barrier(world);
//	MPI_Gatherv(lField, lFSize, MPIField, NULL,
//	 NULL, NULL, NULL, master, world);
//	delete lField;
//}
//
//void
//MPI::setPartitionSizes(
// const vector<vector<size_t> >& partId) {
//	// Sets partition sizes and offsets.
//	pSize = new Math::Int[nTasks];
//	pOffset = new Math::Int [nTasks];
//	fSize = new Math::Int [nTasks];
//	fOffset = new Math::Int [nTasks];
//	for (Math::Int i = 0; i < nTasks; i++) {
//		pSize[i] = partId[i].size();
//		fSize[i] = partId[i].size() * np;
//		if (i == 0) {
//			pOffset[i] = 0;
//			fOffset[i] = 0;
//		} else {
//			pOffset[i] = pSize[i-1] + pOffset[i-1];
//			fOffset[i] = fSize[i-1] + fOffset[i-1];
//		}
//	}
//	setLocalSizeAndOffset(pSize[task], pOffset[task]);
//}
//
//size_t
//MPI::getLocalOffset() const {
//	assert(pOffset != NULL);
//	return pOffset[task];
//}
//
//size_t
//MPI::getLocalSize() const {
//	assert(pSize != NULL);
//	return pSize[task];
//}
//
//void
//MPI::syncNeighbourFields(
// Math::Real* nEx, Math::Real* nEy, Math::Real* nEz,
// Math::Real* nHx, Math::Real* nHy, Math::Real* nHz,
// const Math::Real* Ex, const Math::Real* Ey, const Math::Real* Ez,
// const Math::Real* Hx, const Math::Real* Hy,	const Math::Real* Hz) const {
//	for (Math::Int t = 0; t < nTasks; t++) {
//		Math::Real *sField;
//		sField = new Math::Real[nDofSize[t][task]];
//		// Packs fields for task t.
//		if (t != task) {
//			size_t fInd = 0;
//			for (Math::Int i = 0; i < nNeighId[t]; i++) {
//				size_t id = neighId[t][i];
//				if (isLocalId(id)) {
//					size_t fp = getRelPosOfId(id) * np;
//					for (size_t j = 0; j < np; j++) {
//						sField[fInd++] = Ex[fp];
//						sField[fInd++] = Ey[fp];
//						sField[fInd++] = Ez[fp];
//						sField[fInd++] = Hx[fp];
//						sField[fInd++] = Hy[fp];
//						sField[fInd++] = Hz[fp];
//						fp++;
//					}
//				}
//			}
//		}
//		// Gathers fields.
//		size_t fSize = nNeighId[t] * np;
//		Math::Real *rField;
//		rField = new Math::Real[fSize*6];
//		MPI_Gatherv(sField, nDofSize[t][task], MPI_DOUBLE,
//		 rField, nDofSize[t], nDofOffset[t], MPI_DOUBLE, t, world);
//		// Unpacks fields if this process is the receiver.
//		if (task == t) {
//			size_t nn = (size_t) nNeighId[t];
//			for (size_t i = 0; i < nn; i++) {
//				size_t fInd = i * np * 6;
//				size_t frp = neighElemRP[i] * np;
//				assert((size_t) neighElemRP[i*np] < nn);
//				for (size_t j = 0; j < np; j++) {
//					assert(0 <= frp && frp < nn * np);
//					nEx[frp] = rField[fInd++];
//					nEy[frp] = rField[fInd++];
//					nEz[frp] = rField[fInd++];
//					nHx[frp] = rField[fInd++];
//					nHy[frp] = rField[fInd++];
//					nHz[frp] = rField[fInd++];
//					frp++;
//				}
//			}
//		}
//		// Frees all memory.
//		delete sField;
//		delete rField;
//	}
//}
//
//void
//MPI::initNeighbourFields(const vector<size_t>& nIds) {
//	assert(nNeighId == NULL);
//	assert(neighId == NULL);
//	// Gathers number of neighbours.
//	nNeighId = new Math::Int[nTasks];
//	size_t nNeigh = nIds.size();
//	MPI_Allgather(&nNeigh, 1, MPI_UNSIGNED,
//     nNeighId, 1, MPI_UNSIGNED, world);
//	// Gathers neighbour Ids.
//	// Gathers data.
//	size_t totalNumberOfNeigh = 0;
//	for (Math::Int t = 0; t < nTasks; t++) {
//		totalNumberOfNeigh += nNeighId[t];
//	}
//	Math::Int auxOffset[nTasks];
//	auxOffset[0] = 0;
//	for (Math::Int t = 1; t < nTasks; t++) {
//		auxOffset[t] = nNeighId[t-1] + auxOffset[t-1];
//	}
//	Math::Int neighIdGlobal[totalNumberOfNeigh]; // Aux. receiver buffer.
//	Math::Int neighIdsBuffer[nNeigh];
//	for (size_t i = 0; i < nNeigh; i++) {
//		neighIdsBuffer[i] = nIds[i];
//	}
//	MPI_Allgatherv(neighIdsBuffer, nIds.size(), MPI_UNSIGNED,
//	 neighIdGlobal, nNeighId, auxOffset, MPI_UNSIGNED, world);
//	// Unpacks buffer.
//	neighId = new size_t*[nTasks];
//	for (Math::Int t = 0; t < nTasks; t++) {
//		neighId[t] = new size_t[nNeighId[t]];
//		for (Math::Int i = 0; i < nNeighId[t]; i++) {
//			neighId[t][i] = neighIdGlobal[auxOffset[t] + i];
//		}
//	}
//	// Stores task relational data. Sizes and offsets.
//	assert(neighFSize == NULL);
//	assert(neighFOffset == NULL);
//	neighFSize = new Math::Int *[nTasks];
//	neighFOffset = new Math::Int *[nTasks];
//	nDofSize = new Math::Int *[nTasks];
//	nDofOffset = new Math::Int *[nTasks];
//	for (Math::Int i = 0; i < nTasks; i++) {
//		neighFSize[i] = new Math::Int[nTasks];
//		neighFOffset[i] = new Math::Int[nTasks];
//		nDofSize[i] = new Math::Int[nTasks];
//		nDofOffset[i] = new Math::Int[nTasks];
//	}
//	for (Math::Int t = 0; t < nTasks; t++) {
//		// Counts Ids coming from each task.
//		for (Math::Int i = 0; i < nTasks; i++) {
//			neighFSize[t][i] = 0;
//			nDofSize[t][i] = 0;
//		}
//		for (Math::Int i = 0; i < nNeighId[t]; i++) {
//			Math::Int idTask = getTaskOfId(neighId[t][i]);
//			neighFSize[t][idTask] += np;
//			nDofSize[t][idTask] += np * 6;
//		}
//		// Inits offsets.
//		for (Math::Int i = 0; i < nTasks; i++) {
//			if (i == 0) {
//				neighFOffset[t][0] = 0;
//				nDofOffset[t][0] = 0;
//			} else {
//				neighFOffset[t][i] =
//				 neighFOffset[t][i-1] + neighFSize[t][i-1];
//				nDofOffset[t][i] =
//				 nDofOffset[t][i-1] + nDofSize[t][i-1];
//			}
//		}
//	}
//	// Sends relative positions of neigh nodes in local node.
//	Math::Int **neighSize, **neighOffset;
//	neighSize = new Math::Int*[nTasks];
//	neighOffset = new Math::Int*[nTasks];
//	for (Math::Int t = 0; t < nTasks; t++) {
//		neighSize[t] = new Math::Int[nTasks];
//		neighOffset[t] = new Math::Int[nTasks];
//		// Counts Ids coming from each task.
//		for (Math::Int i = 0; i < nTasks; i++) {
//			neighSize[t][i] = neighFSize[t][i] / np;
//			neighOffset[t][i] = neighFOffset[t][i] / np;
//		}
//	}
//	// Sends relative positions.
//	for (Math::Int t = 0; t < nTasks; t++) {
//		size_t *srp;
//		size_t nFSize = neighFSize[t][task];
//		srp = new size_t[nFSize];
//		if (t != task) {
//			size_t k = 0;
//			for (Math::Int i = 0; i < nNeighId[t]; i++) {
//				size_t id = neighId[t][i];
//				if (isLocalId(id)) {
//					srp[k++] = i;
//				}
//			}
//		} else {
//			neighElemRP = new size_t[nNeighId[t]];
//		}
//		MPI_Gatherv(srp, neighSize[t][task], MPI_UNSIGNED,
//		 neighElemRP, neighSize[t], neighOffset[t], MPI_UNSIGNED,
//		 t, world);
//		if (t != task) {
//			delete srp;
//		}
//	}
//	// Frees memory.
//	delete neighOffset;
//	delete neighSize;
//	// Performs several checkings (only when DBG)
//	if (nTasks == 2) {
//		assert(checkNNeighCoherence(nNeighId));
//	}
//	assert(checkVectorsAreEqual(nNeighId[task], neighId[task], nIds));
//	assert(checkNeighFSizes());
//}
//
//void
//MPI::printInfo() const {
//	cout << " --- CommMPI info --- " << endl;
//	printf ("Number of tasks= %d My rank= %d\n", nTasks, task);
//	cout << " Local Size: " << getLocalSize() << endl;
//	cout << " Global Size: " << getGlobalSize() << endl;
//	if (pSize != NULL) {
//		cout << "pSize: ";
//		for (Math::Int t = 0; t < nTasks; t++) {
//			cout << pSize[t] << " ";
//		}
//		cout << endl;
//	}
//	if (pOffset != NULL) {
//		cout << "pOffset: ";
//		for (Math::Int t = 0; t < nTasks; t++) {
//			cout << pOffset[t] << " ";
//		}
//		cout << endl;
//	}
//	if (nNeighId != NULL) {
//		cout << "nneighId: " << endl;
//		for (Math::Int t = 0; t < nTasks; t++) {
//			cout << "task #" << t << ": " << nNeighId[t] << endl;
//		}
//	}
//	if (neighId != NULL) {
//		cout << "neighId" << endl;
//		for (Math::Int t = 0; t < nTasks; t++) {
//			cout << "task #" << t << ": ";
//			for (Math::Int i = 0; i < nNeighId[t]; i++) {
//				cout << neighId[t][i] << " ";
//			}
//			cout << endl;
//		}
//	}
//	if (neighFSize != NULL) {
//		cout << "neighFSize" << endl;
//		for (Math::Int t = 0; t < nTasks; t++) {
//			for (Math::Int s = 0; s < nTasks; s++) {
//				cout << neighFSize[t][s] << " ";
//			}
//			cout << endl;
//		}
//		cout << "neighFOffset" << endl;
//		for (Math::Int t = 0; t < nTasks; t++) {
//			for (Math::Int s = 0; s < nTasks; s++) {
//				cout << neighFOffset[t][s] << " ";
//			}
//			cout << endl;
//		}
//		cout << "nDofSize" << endl;
//		for (Math::Int t = 0; t < nTasks; t++) {
//			for (Math::Int s = 0; s < nTasks; s++) {
//				cout << nDofSize[t][s] << " ";
//			}
//			cout << endl;
//		}
//		cout << "nDofOffset" << endl;
//		for (Math::Int t = 0; t < nTasks; t++) {
//			for (Math::Int s = 0; s < nTasks; s++) {
//				cout << nDofOffset[t][s] << " ";
//			}
//			cout << endl;
//		}
//	}
//	printOrderingInfo();
//	cout << " --- End of CommMPI info ---" << endl;
//}
//
//void
//MPI::packFields(
// Field *field,
// const Math::Real *Ex, const Math::Real *Ey, const Math::Real *Ez,
// const Math::Real *Hx, const Math::Real *Hy, const Math::Real *Hz,
// const size_t fSize) const {
//	for (size_t i = 0; i < fSize; i++) {
//		field[i].Ex = Ex[i];
//		field[i].Ey = Ey[i];
//		field[i].Ez = Ez[i];
//		field[i].Hx = Hx[i];
//		field[i].Hy = Hy[i];
//		field[i].Hz = Hz[i];
//	}
//}
//
//void
//MPI::unpackFields(
// Math::Real *Ex, Math::Real *Ey, Math::Real *Ez,
// Math::Real *Hx, Math::Real *Hy, Math::Real *Hz,
// const Field *field, const size_t fSize) const {
//	for (size_t i = 0; i < fSize; i++) {
//		Ex[i] = field[i].Ex;
//		Ey[i] = field[i].Ey;
//		Ez[i] = field[i].Ez;
//		Hx[i] = field[i].Hx;
//		Hy[i] = field[i].Hy;
//		Hz[i] = field[i].Hz;
//	}
//}
//
//void
//MPI::commitFieldStruct() {
//	const Math::Int count = 7;
//	Math::Int lengths[count] = {1, 1, 1, 1, 1, 1, 1};
//	Math::Int iTS, dTS;
//	MPI_Type_size(MPI_UNSIGNED, &iTS);
//	MPI_Type_size(MPI_DOUBLE, &dTS);
//	MPI_Aint disp[count] = { 0, iTS, iTS+dTS, iTS+2*dTS,
//	 iTS+3*dTS, iTS+4*dTS, iTS+5*dTS};
//	MPI_Datatype types[count] = {
//	 MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
//	 MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
//	MPI_Type_create_struct(count,lengths,disp,types, &MPIField);
//	MPI_Type_commit(&MPIField);
//}
//
//size_t
//MPI::getLocalFieldSize() const {
//	assert(fSize != NULL);
//	return fSize[task];
//}
//
//Math::Int*
//MPI::getFieldOffsets() const{
//	assert(fOffset != NULL);
//	return fOffset;
//}
//
//Math::Int*
//MPI::getFieldSizes() const{
//	assert(fSize != NULL);
//	return fSize;
//}
//
//
//Math::Int
//MPI::getTaskOfId(const size_t id) const {
//	size_t rp = getGlobalRelPosOfId(id);
//	for (Math::Int t = 0; t < (nTasks - 1); t++) {
//		if (rp < (size_t) pOffset[t+1]) {
//			return t;
//		}
//	}
//	return (nTasks - 1);
//}
//
//bool
//MPI::checkNNeighCoherence(Math::Int* nneigh) const {
//	assert(nneigh != NULL);
//	return (nneigh[0] == nneigh[1]);
//}
//
//bool
//MPI::checkVectorsAreEqual(
// const size_t vSize,
// const size_t* v,
// const vector<size_t>& nIds) const {
//	assert(vSize == nIds.size());
//	bool ok = true;
//	for (size_t i = 0; i < vSize; i++) {
//		if (v[i] != nIds[i]) {
//			cout << v[i] << " is not equal to " << nIds[i] << endl;
//			ok = false;
//		}
//	}
//	if (!ok) {
//		cout << "ERROR@checkVectorsAreEqual @ Task: " << task << endl;
//		cout << "Above values of gathered neighbour ids are not equal"
//		 << endl;
//	}
//	return ok;
//}
//
//Math::Real
//MPI::reduceToGlobalMinimum(Math::Real val) const {
//	Math::Real res;
//	static const Math::Int count = 1;
//	MPI_Allreduce(&val, &res, count,
//	 MPI_DOUBLE, MPI_MIN, world);
//	return res;
//}
//
//Math::Int
//MPI::countTasksInLocalHost() const {
//	// Counts number of tasks running on this host.
//	char localHostName[MPI_MAX_PROCESSOR_NAME];
//	Math::Int localHostNameLen;
//	MPI_Get_processor_name(localHostName, &localHostNameLen);
//	char hostName[nTasks][MPI_MAX_PROCESSOR_NAME];
//	Math::Int hostNameLen[nTasks];
//	MPI_Allgather(
//	 &localHostNameLen, 1, MPI_INT, hostNameLen, 1, MPI_INT,
//		world);
//	MPI_Allgather(
//	 localHostName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, hostName,
//	 MPI_MAX_PROCESSOR_NAME, MPI_CHAR, world);
//	Math::Int res = 0;
//	for (Math::Int t = 0; t < nTasks; t++) {
//		bool isLocalHost = true;
//		if (localHostNameLen != hostNameLen[t]) {
//			isLocalHost &= false;
//		}
//		for (Math::Int i = 0; i < localHostNameLen; i++) {
//			isLocalHost &= (localHostName[i] == hostName[t][i]);
//		}
//		if (isLocalHost) {
//			res++;
//		}
//	}
//	assert(res > 0);
//	return res;
//}
//
//Math::Int
//MPI::getNumOfTasksOnThisHost() const {
//	return nTasksInHost;
//}
//
//bool
//MPI::checkNeighFSizes() const {
//	bool sizesOk = true;
//	for (Math::Int t = 0; t < nTasks; t++) {
//		size_t neighFSum = 0;
//		for (Math::Int s = 0; s < nTasks; s++) {
//			neighFSum += neighFSize[t][s];
//		}
//		sizesOk &= (neighFSum == nNeighId[t] * np);
//	}
//	if (!sizesOk) {
//		cerr << endl << "neighFSizes are inconsistent with nNeighId" << endl;
//	}
//	bool diagOk = true;
//	for (Math::Int t = 0; t < nTasks; t++) {
//			diagOk &= (neighFSize[t][t] == 0);
//	}
//	if (!diagOk) {
//		printInfo();
//		cerr << endl << "neighFSizes has non-zero diag." << endl;
//		cerr << endl << "A task is neighbouring itself." << endl;
//	}
//	bool symmetryOk = true;
//	DynMatrix<Math::Int> aux(nTasks, nTasks, neighFSize);
//	symmetryOk = aux.isSymmetric();
//	if (!symmetryOk) {
//		printInfo();
//		cerr << endl << "ERROR @ checkNeighSizes()" << endl;
//		cerr << endl << "neighFSizes is not symmetric." << endl;
//	}
//	return (sizesOk && diagOk && symmetryOk);
//}
//
//vector<size_t>
//MPI::getThreadsOfTasks() const {
//	size_t *taskThreads;
//	taskThreads = new size_t[nTasks];
//#ifdef USE_OPENMP
//	size_t localThreads = omp_get_max_threads();
//#else
//	size_t localThreads = 1;
//#endif
//	MPI_Allgather(&localThreads, 1, MPI_UNSIGNED,
//     taskThreads, 1, MPI_UNSIGNED, world);
//	vector<size_t> res;
//	res.reserve(nTasks);
//	for (Math::Int t = 0; t < nTasks; t++) {
//		res.push_back(taskThreads[t]);
//	}
//	delete taskThreads;
//	return res;
//}
//
//}
//}
