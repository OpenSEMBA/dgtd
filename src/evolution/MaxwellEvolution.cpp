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

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << std::endl;
#endif

		if (pd_.model.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != pd_.model.getTotalFieldScatteredFieldToMarker().end()) {
			srcmngr.initTFSFPreReqs(pd_.model.getConstMesh(), pd_.model.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));
			auto globalTFSFfes{ srcmngr.getGlobalTFSFSpace() };			
			Model tfsfModel = Model(*globalTFSFfes->GetMesh(), GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
			ProblemDescription tfsfPD(tfsfModel, pd_.probes, pd_.sources, pd_.opts);
			DGOperatorFactory tfsfFactory(tfsfPD, *globalTFSFfes);

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Assembling TFSF Inverse Mass Operators" << std::endl;
#endif

			MInvTFSF_ = tfsfFactory.buildMaxwellInverseMassMatrixOperator();

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass Zero-Normal Operators" << std::endl;
#endif

			MP_GTFSF_ = tfsfFactory.buildMaxwellZeroNormalOperator();

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass One-Normal Operators" << std::endl;
#endif
			MFN_GTFSF_ = tfsfFactory.buildMaxwellOneNormalOperator();

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass Two-Normal Operators" << std::endl;
#endif

			MFNN_GTFSF_ = tfsfFactory.buildMaxwellTwoNormalOperator();

		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling Standard Inverse Mass Operators" << std::endl;
#endif

		DGOperatorFactory dgFactory(pd_, fes);

		if (pd_.model.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
#endif
			MPB_ = dgFactory.buildMaxwellIntBdrZeroNormalOperator();

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
#endif
			MFNB_ = dgFactory.buildMaxwellIntBdrOneNormalOperator();

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
#endif

			MFNNB_ = dgFactory.buildMaxwellIntBdrTwoNormalOperator();

		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
#endif

		MS_ = dgFactory.buildMaxwellDirectionalOperator();

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
#endif
		
		MP_ = dgFactory.buildMaxwellZeroNormalOperator();

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
#endif

		MFN_ = dgFactory.buildMaxwellOneNormalOperator();

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
#endif

		MFNN_ = dgFactory.buildMaxwellTwoNormalOperator();

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Operator assembly finished" << std::endl;
		std::cout << std::endl;
#endif

	}

void MaxwellEvolution::Mult(const Vector& in, Vector& out) const
{

	std::array<Vector, 3> eOld, hOld;
	std::array<ParGridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetDataAndSize(in.GetData() + d * fes_.GetNDofs(), fes_.GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d + 3) * fes_.GetNDofs(), fes_.GetNDofs());
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

		MFN_[H][E][y]->AddMult(eOld[z], hNew[x], 1.0);
		MFN_[H][E][z]->AddMult(eOld[y], hNew[x], -1.0);
		MFN_[E][H][y]->AddMult(hOld[z], eNew[x], -1.0);
		MFN_[E][H][z]->AddMult(hOld[y], eNew[x], 1.0);

		//Upwind

		MFNN_[H][H][X][x]->AddMult(hOld[X], hNew[x], 1.0);
		MFNN_[H][H][Y][x]->AddMult(hOld[Y], hNew[x], 1.0);
		MFNN_[H][H][Z][x]->AddMult(hOld[Z], hNew[x], 1.0);
		MP_[H]->AddMult(hOld[x], hNew[x], -1.0);

		MFNN_[E][E][Y][x]->AddMult(eOld[Y], eNew[x], 1.0);
		MFNN_[E][E][X][x]->AddMult(eOld[X], eNew[x], 1.0);
		MFNN_[E][E][Z][x]->AddMult(eOld[Z], eNew[x], 1.0);
		MP_[E]->AddMult(eOld[x], eNew[x], -1.0);
		

		if (pd_.model.getInteriorBoundaryToMarker().size()) {

			MFNB_[H][E][y]->AddMult(eOld[z], hNew[x]);
			MFNB_[H][E][z]->AddMult(eOld[y], hNew[x], -1.0);
			MFNB_[E][H][y]->AddMult(hOld[z], eNew[x], -1.0);
			MFNB_[E][H][z]->AddMult(hOld[y], eNew[x]);

			//Upwind

			MFNNB_[H][H][X][x]->AddMult(hOld[X], hNew[x]);
			MFNNB_[H][H][Y][x]->AddMult(hOld[Y], hNew[x]);
			MFNNB_[H][H][Z][x]->AddMult(hOld[Z], hNew[x]);
			MPB_[H]->AddMult(hOld[x], hNew[x], -1.0);

			MFNNB_[E][E][X][x]->AddMult(eOld[X], eNew[x]);
			MFNNB_[E][E][Y][x]->AddMult(eOld[Y], eNew[x]);
			MFNNB_[E][E][Z][x]->AddMult(eOld[Z], eNew[x]);
			MPB_[E]->AddMult(eOld[x], eNew[x], -1.0);
			
		}
	}

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<TotalField*>(source.get())) {

			auto func{ evalTimeVarFunction(GetTime(),srcmngr_) };

			std::array<ParGridFunction, 3> eTemp, hTemp;

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
}

}

