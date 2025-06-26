#include "MaxwellEvolution.h"

#include <chrono>
#include <algorithm>
#include <iostream>
#include <string>

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

void evalConductivity(const Vector& cond, const Vector& in, Vector& out)
{
	for (auto v{ 0 }; v < cond.Size(); v++) {
		out[v] -= cond[v] * in[v];
	}
}

void changeSignOfFieldGridFuncs(FieldGridFuncs& gfs)
{
	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			gfs[f][d] *= -1.0;
		}
	}
}

static void loadFluxVector(const ParGridFunction& local, const Vector& nbr, Vector& flux)
{
	flux.SetSize(local.Size() + nbr.Size());
	for (auto v = 0; v < local.Size(); v++){
		flux[v] = local[v];
	}
	for (auto v = 0; v < nbr.Size(); v++){
		flux[v+local.Size()] = nbr[v];
	}
}

MaxwellEvolution::MaxwellEvolution(
	ProblemDescription& pd, ParFiniteElementSpace& fes, SourcesManager& srcmngr) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	pd_(pd),
	fes_(fes),
	srcmngr_(srcmngr)
{
#ifdef SHOW_TIMER_INFORMATION
		auto startTime{ std::chrono::high_resolution_clock::now() };
#endif
		fes_.ExchangeFaceNbrData();

#ifdef SHOW_TIMER_INFORMATION
	if (Mpi::WorldRank() == 0){
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << std::endl;
	}
#endif

		if (pd_.model.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != pd_.model.getTotalFieldScatteredFieldToMarker().end()) {
			srcmngr.initTFSFPreReqs(pd_.model.getConstMesh(), pd_.model.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));
			auto globalTFSFfes{ srcmngr.getGlobalTFSFSpace() };			
			Model tfsfModel(*globalTFSFfes->GetMesh(), GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
			ProblemDescription tfsfPD(tfsfModel, pd_.probes, pd_.sources, pd_.opts);
			DGOperatorFactory<FiniteElementSpace> tfsfFactory(tfsfPD, *globalTFSFfes);

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Assembling TFSF Inverse Mass Operators" << std::endl;
		}
#endif

			MInvTFSF_ = tfsfFactory.buildMaxwellInverseMassMatrixOperator<BilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass Zero-Normal Operators" << std::endl;
		}
#endif

			MP_GTFSF_ = tfsfFactory.buildMaxwellZeroNormalOperator<BilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass One-Normal Operators" << std::endl;
		}
#endif
			MFN_GTFSF_ = tfsfFactory.buildMaxwellOneNormalOperator<BilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass Two-Normal Operators" << std::endl;
		}
#endif

			MFNN_GTFSF_ = tfsfFactory.buildMaxwellTwoNormalOperator<BilinearForm>();

		}

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Operators" << std::endl;
		}
#endif

		DGOperatorFactory<ParFiniteElementSpace> dgFactory(pd_, fes_);

		if (pd_.model.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
		}
#endif
			MPB_ = dgFactory.buildMaxwellIntBdrZeroNormalOperator<ParBilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
		}
#endif
			MFNB_ = dgFactory.buildMaxwellIntBdrOneNormalOperator<ParBilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
		}
#endif

			MFNNB_ = dgFactory.buildMaxwellIntBdrTwoNormalOperator<ParBilinearForm>();

		}

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
		}
#endif

		MS_ = dgFactory.buildMaxwellDirectionalOperator<ParBilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
		}
#endif
		
		MP_ = dgFactory.buildMaxwellZeroNormalOperator<ParBilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
		}
#endif

		MFN_ = dgFactory.buildMaxwellOneNormalOperator<ParBilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
	}
#endif

		MFNN_ = dgFactory.buildMaxwellTwoNormalOperator<ParBilinearForm>();

#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Operator assembly finished" << std::endl;
			std::cout << std::endl;
		}
#endif

}

void MaxwellEvolution::Mult(const Vector& in, Vector& out) const
{

	std::array<ParGridFunction, 3> eOld, hOld, eNew, hNew, eFlux, hFlux;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetSpace(&fes_);
		hOld[d].SetSpace(&fes_);
		eOld[d].SetDataAndSize(in.GetData() + d * fes_.GetNDofs(), fes_.GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d + 3) * fes_.GetNDofs(), fes_.GetNDofs());
		eOld[d].ExchangeFaceNbrData();
		hOld[d].ExchangeFaceNbrData();
		auto faceNbrField = eOld[d].FaceNbrData();
		loadFluxVector(eOld[d], faceNbrField, eFlux[d]);
		faceNbrField = hOld[d].FaceNbrData();
		loadFluxVector(hOld[d], faceNbrField, hFlux[d]);

		eNew[d].SetSpace(&fes_);
		hNew[d].SetSpace(&fes_);
		eNew[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hNew[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
		eNew[d] = 0.0;
		hNew[d] = 0.0;
	}

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		//Centered

		MS_[H][y]->AddMult(eOld[z], hNew[x], -1.0);
		MS_[H][z]->AddMult(eOld[y], hNew[x]);
		MS_[E][y]->AddMult(hOld[z], eNew[x]);
		MS_[E][z]->AddMult(hOld[y], eNew[x], -1.0);

		MFN_[H][E][y]->AddMult(eFlux[z], hNew[x], 1.0);
		MFN_[H][E][z]->AddMult(eFlux[y], hNew[x], -1.0);
		MFN_[E][H][y]->AddMult(hFlux[z], eNew[x], -1.0);
		MFN_[E][H][z]->AddMult(hFlux[y], eNew[x], 1.0);

		//Upwind

		MFNN_[H][H][X][x]->AddMult(hFlux[X], hNew[x], 1.0);
		MFNN_[H][H][Y][x]->AddMult(hFlux[Y], hNew[x], 1.0);
		MFNN_[H][H][Z][x]->AddMult(hFlux[Z], hNew[x], 1.0);
		MP_[H]->AddMult(hFlux[x], hNew[x], -1.0);

		MFNN_[E][E][Y][x]->AddMult(eFlux[Y], eNew[x], 1.0);
		MFNN_[E][E][X][x]->AddMult(eFlux[X], eNew[x], 1.0);
		MFNN_[E][E][Z][x]->AddMult(eFlux[Z], eNew[x], 1.0);
		MP_[E]->AddMult(eFlux[x], eNew[x], -1.0);
		

		if (pd_.model.getInteriorBoundaryToMarker().size()) {

			MFNB_[H][E][y]->AddMult(eFlux[z], hNew[x]);
			MFNB_[H][E][z]->AddMult(eFlux[y], hNew[x], -1.0);
			MFNB_[E][H][y]->AddMult(hFlux[z], eNew[x], -1.0);
			MFNB_[E][H][z]->AddMult(hFlux[y], eNew[x]);

			//Upwind

			MFNNB_[H][H][X][x]->AddMult(hFlux[X], hNew[x]);
			MFNNB_[H][H][Y][x]->AddMult(hFlux[Y], hNew[x]);
			MFNNB_[H][H][Z][x]->AddMult(hFlux[Z], hNew[x]);
			MPB_[H]->AddMult(hFlux[x], hNew[x], -1.0);

			MFNNB_[E][E][X][x]->AddMult(eFlux[X], eNew[x]);
			MFNNB_[E][E][Y][x]->AddMult(eFlux[Y], eNew[x]);
			MFNNB_[E][E][Z][x]->AddMult(eFlux[Z], eNew[x]);
			MPB_[E]->AddMult(eFlux[x], eNew[x], -1.0);
			
		}
	}

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<TotalField*>(source.get()) && srcmngr_.getGlobalTFSFSpace() != nullptr) {

			auto func{ evalTimeVarFunction(GetTime(),srcmngr_) };

			std::array<GridFunction, 3> eTemp, hTemp;

			for (int d = X; d <= Z; d++) {
				eTemp[d].SetSpace(srcmngr_.getGlobalTFSFSpace());
				hTemp[d].SetSpace(srcmngr_.getGlobalTFSFSpace());
				eTemp[d] = 0.0;
				hTemp[d] = 0.0;
			}

			for (int x = X; x <= Z; x++) {
				int y = (x + 1) % 3;
				int z = (x + 2) % 3;

				MaxwellTransferMap eMap(eTemp[x], eNew[x]);
				MaxwellTransferMap hMap(hTemp[x], hNew[x]);

				//Centered

				MFN_GTFSF_[H][E][y]->Mult(func[E][z], hTemp[x]);
				eMap.TransferSub(hTemp[x], hNew[x]);
				MFN_GTFSF_[H][E][z]->Mult(func[E][y], hTemp[x]);
				eMap.TransferAdd(hTemp[x], hNew[x]);
				MFN_GTFSF_[E][H][y]->Mult(func[H][z], eTemp[x]);
				eMap.TransferAdd(eTemp[x], eNew[x]);
				MFN_GTFSF_[E][H][z]->Mult(func[H][y], eTemp[x]);
				eMap.TransferSub(eTemp[x], eNew[x]);

				//Upwind

				MFNN_GTFSF_[H][H][X][x]->Mult(func[H][X], hTemp[x]);
				hMap.TransferSub(hTemp[x], hNew[x]);
				MFNN_GTFSF_[H][H][Y][x]->Mult(func[H][Y], hTemp[x]);
				hMap.TransferSub(hTemp[x], hNew[x]);
				MFNN_GTFSF_[H][H][Z][x]->Mult(func[H][Z], hTemp[x]);
				hMap.TransferSub(hTemp[x], hNew[x]);
				MP_GTFSF_[H]->Mult(func[H][x], hTemp[x]);
				hMap.TransferAdd(hTemp[x], hNew[x]);

				MFNN_GTFSF_[E][E][X][x]->Mult(func[E][X], eTemp[x]);
				eMap.TransferSub(eTemp[x], eNew[x]);
				MFNN_GTFSF_[E][E][Y][x]->Mult(func[E][Y], eTemp[x]);
				eMap.TransferSub(eTemp[x], eNew[x]);
				MFNN_GTFSF_[E][E][Z][x]->Mult(func[E][Z], eTemp[x]);
				eMap.TransferSub(eTemp[x], eNew[x]);
				MP_GTFSF_[E]->Mult(func[E][x], eTemp[x]);
				eMap.TransferAdd(eTemp[x], eNew[x]);

			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

}

