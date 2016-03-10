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
///*
// * CommunicationsMPI.cpp
// *
// *  Created on: Apr 17, 2013
// *      Author: luis
// */
//
//#include "CommMPI.h"
//
//CommMPI::CommMPI() {
//	world = MPI_COMM_WORLD;
//	Int rc = MPI_Init(NULL, NULL);
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
//CommMPI::~CommMPI() {
//	MPI_Finalize();
//}
//
//Int
//CommMPI::getNumberOfTasks() const {
//	return nTasks;
//}
//
//Int
//CommMPI::getTask() const {
//	return task;
//}
//
//bool
//CommMPI::isMaster() const {
//	if (task == master) {
//		return true;
//	}
//	return false;
//}
//
//void
//CommMPI::gatherFieldsMaster(
// Real* Ex, Real* Ey, Real* Ez,
// Real* Hx, Real* Hy, Real* Hz,
// const Real* lEx, const Real* lEy, const Real* lEz,
// const Real* lHx, const Real* lHy, const Real* lHz) const {
//	assert(isMaster());
//	// Packs fields in field struct.
//	UInt lFieldSize = getLocalFieldSize();
//	Field_s *lField;
//	lField = new Field_s[lFieldSize];
//	packFields(lField, lEx,lEy,lEz,lHx,lHy,lHz, lFieldSize);
//	// Allocates space for general field.
//	UInt gFieldSize = getGlobalSize() * np;
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
//CommMPI::gatherFieldsSlave(
// const Real* lEx, const Real* lEy, const Real* lEz,
// const Real* lHx, const Real* lHy, const Real* lHz) const {
//	assert(!isMaster());
//	// Packs fields in field struct.
//	UInt lFSize = getLocalFieldSize();
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
//CommMPI::setPartitionSizes(
// const vector<vector<UInt> >& partId) {
//	// Sets partition sizes and offsets.
//	pSize = new Int[nTasks];
//	pOffset = new Int [nTasks];
//	fSize = new Int [nTasks];
//	fOffset = new Int [nTasks];
//	for (Int i = 0; i < nTasks; i++) {
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
//UInt
//CommMPI::getLocalOffset() const {
//	assert(pOffset != NULL);
//	return pOffset[task];
//}
//
//UInt
//CommMPI::getLocalSize() const {
//	assert(pSize != NULL);
//	return pSize[task];
//}
//
//void
//CommMPI::syncNeighbourFields(
// Real* nEx, Real* nEy, Real* nEz,
// Real* nHx, Real* nHy, Real* nHz,
// const Real* Ex, const Real* Ey, const Real* Ez,
// const Real* Hx, const Real* Hy,	const Real* Hz) const {
//	for (Int t = 0; t < nTasks; t++) {
//		Real *sField;
//		sField = new Real[nDofSize[t][task]];
//		// Packs fields for task t.
//		if (t != task) {
//			UInt fInd = 0;
//			for (Int i = 0; i < nNeighId[t]; i++) {
//				UInt id = neighId[t][i];
//				if (isLocalId(id)) {
//					UInt fp = getRelPosOfId(id) * np;
//					for (UInt j = 0; j < np; j++) {
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
//		UInt fSize = nNeighId[t] * np;
//		Real *rField;
//		rField = new Real[fSize*6];
//		MPI_Gatherv(sField, nDofSize[t][task], MPI_DOUBLE,
//		 rField, nDofSize[t], nDofOffset[t], MPI_DOUBLE, t, world);
//		// Unpacks fields if this process is the receiver.
//		if (task == t) {
//			UInt nn = (UInt) nNeighId[t];
//			for (UInt i = 0; i < nn; i++) {
//				UInt fInd = i * np * 6;
//				UInt frp = neighElemRP[i] * np;
//				assert((UInt) neighElemRP[i*np] < nn);
//				for (UInt j = 0; j < np; j++) {
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
//CommMPI::initNeighbourFields(const vector<UInt>& nIds) {
//	assert(nNeighId == NULL);
//	assert(neighId == NULL);
//	// Gathers number of neighbours.
//	nNeighId = new Int[nTasks];
//	UInt nNeigh = nIds.size();
//	MPI_Allgather(&nNeigh, 1, MPI_UNSIGNED,
//     nNeighId, 1, MPI_UNSIGNED, world);
//	// Gathers neighbour Ids.
//	// Gathers data.
//	UInt totalNumberOfNeigh = 0;
//	for (Int t = 0; t < nTasks; t++) {
//		totalNumberOfNeigh += nNeighId[t];
//	}
//	Int auxOffset[nTasks];
//	auxOffset[0] = 0;
//	for (Int t = 1; t < nTasks; t++) {
//		auxOffset[t] = nNeighId[t-1] + auxOffset[t-1];
//	}
//	Int neighIdGlobal[totalNumberOfNeigh]; // Aux. receiver buffer.
//	Int neighIdsBuffer[nNeigh];
//	for (UInt i = 0; i < nNeigh; i++) {
//		neighIdsBuffer[i] = nIds[i];
//	}
//	MPI_Allgatherv(neighIdsBuffer, nIds.size(), MPI_UNSIGNED,
//	 neighIdGlobal, nNeighId, auxOffset, MPI_UNSIGNED, world);
//	// Unpacks buffer.
//	neighId = new UInt*[nTasks];
//	for (Int t = 0; t < nTasks; t++) {
//		neighId[t] = new UInt[nNeighId[t]];
//		for (Int i = 0; i < nNeighId[t]; i++) {
//			neighId[t][i] = neighIdGlobal[auxOffset[t] + i];
//		}
//	}
//	// Stores task relational data. Sizes and offsets.
//	assert(neighFSize == NULL);
//	assert(neighFOffset == NULL);
//	neighFSize = new Int *[nTasks];
//	neighFOffset = new Int *[nTasks];
//	nDofSize = new Int *[nTasks];
//	nDofOffset = new Int *[nTasks];
//	for (Int i = 0; i < nTasks; i++) {
//		neighFSize[i] = new Int[nTasks];
//		neighFOffset[i] = new Int[nTasks];
//		nDofSize[i] = new Int[nTasks];
//		nDofOffset[i] = new Int[nTasks];
//	}
//	for (Int t = 0; t < nTasks; t++) {
//		// Counts Ids coming from each task.
//		for (Int i = 0; i < nTasks; i++) {
//			neighFSize[t][i] = 0;
//			nDofSize[t][i] = 0;
//		}
//		for (Int i = 0; i < nNeighId[t]; i++) {
//			Int idTask = getTaskOfId(neighId[t][i]);
//			neighFSize[t][idTask] += np;
//			nDofSize[t][idTask] += np * 6;
//		}
//		// Inits offsets.
//		for (Int i = 0; i < nTasks; i++) {
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
//	Int **neighSize, **neighOffset;
//	neighSize = new Int*[nTasks];
//	neighOffset = new Int*[nTasks];
//	for (Int t = 0; t < nTasks; t++) {
//		neighSize[t] = new Int[nTasks];
//		neighOffset[t] = new Int[nTasks];
//		// Counts Ids coming from each task.
//		for (Int i = 0; i < nTasks; i++) {
//			neighSize[t][i] = neighFSize[t][i] / np;
//			neighOffset[t][i] = neighFOffset[t][i] / np;
//		}
//	}
//	// Sends relative positions.
//	for (Int t = 0; t < nTasks; t++) {
//		UInt *srp;
//		UInt nFSize = neighFSize[t][task];
//		srp = new UInt[nFSize];
//		if (t != task) {
//			UInt k = 0;
//			for (Int i = 0; i < nNeighId[t]; i++) {
//				UInt id = neighId[t][i];
//				if (isLocalId(id)) {
//					srp[k++] = i;
//				}
//			}
//		} else {
//			neighElemRP = new UInt[nNeighId[t]];
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
//CommMPI::printInfo() const {
//	cout << " --- CommMPI info --- " << endl;
//	printf ("Number of tasks= %d My rank= %d\n", nTasks, task);
//	cout << " Local Size: " << getLocalSize() << endl;
//	cout << " Global Size: " << getGlobalSize() << endl;
//	if (pSize != NULL) {
//		cout << "pSize: ";
//		for (Int t = 0; t < nTasks; t++) {
//			cout << pSize[t] << " ";
//		}
//		cout << endl;
//	}
//	if (pOffset != NULL) {
//		cout << "pOffset: ";
//		for (Int t = 0; t < nTasks; t++) {
//			cout << pOffset[t] << " ";
//		}
//		cout << endl;
//	}
//	if (nNeighId != NULL) {
//		cout << "nneighId: " << endl;
//		for (Int t = 0; t < nTasks; t++) {
//			cout << "task #" << t << ": " << nNeighId[t] << endl;
//		}
//	}
//	if (neighId != NULL) {
//		cout << "neighId" << endl;
//		for (Int t = 0; t < nTasks; t++) {
//			cout << "task #" << t << ": ";
//			for (Int i = 0; i < nNeighId[t]; i++) {
//				cout << neighId[t][i] << " ";
//			}
//			cout << endl;
//		}
//	}
//	if (neighFSize != NULL) {
//		cout << "neighFSize" << endl;
//		for (Int t = 0; t < nTasks; t++) {
//			for (Int s = 0; s < nTasks; s++) {
//				cout << neighFSize[t][s] << " ";
//			}
//			cout << endl;
//		}
//		cout << "neighFOffset" << endl;
//		for (Int t = 0; t < nTasks; t++) {
//			for (Int s = 0; s < nTasks; s++) {
//				cout << neighFOffset[t][s] << " ";
//			}
//			cout << endl;
//		}
//		cout << "nDofSize" << endl;
//		for (Int t = 0; t < nTasks; t++) {
//			for (Int s = 0; s < nTasks; s++) {
//				cout << nDofSize[t][s] << " ";
//			}
//			cout << endl;
//		}
//		cout << "nDofOffset" << endl;
//		for (Int t = 0; t < nTasks; t++) {
//			for (Int s = 0; s < nTasks; s++) {
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
//CommMPI::packFields(
// Field *field,
// const Real *Ex, const Real *Ey, const Real *Ez,
// const Real *Hx, const Real *Hy, const Real *Hz,
// const UInt fSize) const {
//	for (UInt i = 0; i < fSize; i++) {
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
//CommMPI::unpackFields(
// Real *Ex, Real *Ey, Real *Ez,
// Real *Hx, Real *Hy, Real *Hz,
// const Field *field, const UInt fSize) const {
//	for (UInt i = 0; i < fSize; i++) {
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
//CommMPI::commitFieldStruct() {
//	const Int count = 7;
//	Int lengths[count] = {1, 1, 1, 1, 1, 1, 1};
//	Int iTS, dTS;
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
//UInt
//CommMPI::getLocalFieldSize() const {
//	assert(fSize != NULL);
//	return fSize[task];
//}
//
//Int*
//CommMPI::getFieldOffsets() const{
//	assert(fOffset != NULL);
//	return fOffset;
//}
//
//Int*
//CommMPI::getFieldSizes() const{
//	assert(fSize != NULL);
//	return fSize;
//}
//
//
//Int
//CommMPI::getTaskOfId(const UInt id) const {
//	UInt rp = getGlobalRelPosOfId(id);
//	for (Int t = 0; t < (nTasks - 1); t++) {
//		if (rp < (UInt) pOffset[t+1]) {
//			return t;
//		}
//	}
//	return (nTasks - 1);
//}
//
//bool
//CommMPI::checkNNeighCoherence(Int* nneigh) const {
//	assert(nneigh != NULL);
//	return (nneigh[0] == nneigh[1]);
//}
//
//bool
//CommMPI::checkVectorsAreEqual(
// const UInt vSize,
// const UInt* v,
// const vector<UInt>& nIds) const {
//	assert(vSize == nIds.size());
//	bool ok = true;
//	for (UInt i = 0; i < vSize; i++) {
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
//Real
//CommMPI::reduceToGlobalMinimum(Real val) const {
//	Real res;
//	static const Int count = 1;
//	MPI_Allreduce(&val, &res, count,
//	 MPI_DOUBLE, MPI_MIN, world);
//	return res;
//}
//
//Int
//CommMPI::countTasksInLocalHost() const {
//	// Counts number of tasks running on this host.
//	char localHostName[MPI_MAX_PROCESSOR_NAME];
//	Int localHostNameLen;
//	MPI_Get_processor_name(localHostName, &localHostNameLen);
//	char hostName[nTasks][MPI_MAX_PROCESSOR_NAME];
//	Int hostNameLen[nTasks];
//	MPI_Allgather(
//	 &localHostNameLen, 1, MPI_INT, hostNameLen, 1, MPI_INT,
//		world);
//	MPI_Allgather(
//	 localHostName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, hostName,
//	 MPI_MAX_PROCESSOR_NAME, MPI_CHAR, world);
//	Int res = 0;
//	for (Int t = 0; t < nTasks; t++) {
//		bool isLocalHost = true;
//		if (localHostNameLen != hostNameLen[t]) {
//			isLocalHost &= false;
//		}
//		for (Int i = 0; i < localHostNameLen; i++) {
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
//Int
//CommMPI::getNumOfTasksOnThisHost() const {
//	return nTasksInHost;
//}
//
//bool
//CommMPI::checkNeighFSizes() const {
//	bool sizesOk = true;
//	for (Int t = 0; t < nTasks; t++) {
//		UInt neighFSum = 0;
//		for (Int s = 0; s < nTasks; s++) {
//			neighFSum += neighFSize[t][s];
//		}
//		sizesOk &= (neighFSum == nNeighId[t] * np);
//	}
//	if (!sizesOk) {
//		cerr << endl << "neighFSizes are inconsistent with nNeighId" << endl;
//	}
//	bool diagOk = true;
//	for (Int t = 0; t < nTasks; t++) {
//			diagOk &= (neighFSize[t][t] == 0);
//	}
//	if (!diagOk) {
//		printInfo();
//		cerr << endl << "neighFSizes has non-zero diag." << endl;
//		cerr << endl << "A task is neighbouring itself." << endl;
//	}
//	bool symmetryOk = true;
//	DynMatrix<Int> aux(nTasks, nTasks, neighFSize);
//	symmetryOk = aux.isSymmetric();
//	if (!symmetryOk) {
//		printInfo();
//		cerr << endl << "ERROR @ checkNeighSizes()" << endl;
//		cerr << endl << "neighFSizes is not symmetric." << endl;
//	}
//	return (sizesOk && diagOk && symmetryOk);
//}
//
//vector<UInt>
//CommMPI::getThreadsOfTasks() const {
//	UInt *taskThreads;
//	taskThreads = new UInt[nTasks];
//#ifdef USE_OPENMP
//	UInt localThreads = omp_get_max_threads();
//#else
//	UInt localThreads = 1;
//#endif
//	MPI_Allgather(&localThreads, 1, MPI_UNSIGNED,
//     taskThreads, 1, MPI_UNSIGNED, world);
//	vector<UInt> res;
//	res.reserve(nTasks);
//	for (Int t = 0; t < nTasks; t++) {
//		res.push_back(taskThreads[t]);
//	}
//	delete taskThreads;
//	return res;
//}
