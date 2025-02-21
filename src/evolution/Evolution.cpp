#include "Evolution.h"

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

const FieldGridFuncs evalTimeVarFunction(const Time time, SourcesManager& sm)
{
	auto res{ sm.evalTimeVarField(time, sm.getGlobalTFSFSpace()) };
	auto func_g_sf = res;
	sm.markDoFSforTForSF(res, true);
	{
		if (sm.getTFSFSubMesher().getSFSubMesh() != NULL) {
			sm.markDoFSforTForSF(func_g_sf, false);
			for (int f : {E, H}) {
				for (int x{ 0 }; x <= Z; x++) {
					res[f][x] -= func_g_sf[f][x];
					res[f][x] *= 0.5;
				}
			}
		}
	}
	return res;
}

MaxwellEvolution::MaxwellEvolution(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options },
	TFSFOperator_{ 0 },
	globalOperator_{ numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs() }
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

	if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {
		srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));
		auto globalTFSFfes{ srcmngr_.getGlobalTFSFSpace() };
		Model modelGlobal = Model(*globalTFSFfes->GetMesh(), GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Assembling TFSF Inverse Mass Operators" << std::endl;
#endif
		std::array<FiniteElementOperator, 2> MInvTFSF;
		for (auto f : { E, H }) {
			MInvTFSF[f] = buildInverseMassMatrix(f, modelGlobal, *globalTFSFfes);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono:: >
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass Zero-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			MP_GTFSF_[f] = buildByMult(*MInvTFSF[f], *buildZeroNormalOperator(f, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass One-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					MFN_GTFSF_[f][f2][d] = buildByMult(*MInvTFSF[f], *buildOneNormalOperator(f2, {d}, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
				}
			}
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass Two-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						MFNN_GTFSF_[f][f2][d][d2] = buildByMult(*MInvTFSF[f], *buildTwoNormalOperator(f2, {d, d2}, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
					}
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Operators" << std::endl;
#endif

	GlobalIndices globalId(fes_.GetNDofs());
	std::array<FiniteElementOperator, 2> MInv;
	MInv[E] = buildInverseMassMatrix(E, model_, fes_);
	MInv[H] = buildInverseMassMatrix(H, model_, fes_);

	DenseMatrix* denseMat;

	if (model_.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto intBdrFluxOneNormal = buildByMult(*MInv[f], *buildOneNormalIBFIOperator(altField(f), { x }, model_, fes_, opts_), fes_);
				denseMat = intBdrFluxOneNormal->SpMat().ToDenseMatrix();
				*denseMat *= -1.0 + double(f) * 2.0; //anticyclic permutations are negative for E, positive for H.
				globalOperator_.AddSubMatrix(globalId.index[f][y], globalId.index[altField(f)][z], *denseMat);
				delete denseMat;
			}
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto intBdrFluxOneNormal = buildByMult(*MInv[f], *buildOneNormalIBFIOperator(altField(f), { x }, model_, fes_, opts_), fes_);
				denseMat = intBdrFluxOneNormal->SpMat().ToDenseMatrix();
				*denseMat *= 1.0 - double(f) * 2.0; //cyclic permutations are positive for E, negative for H.
				globalOperator_.AddSubMatrix(globalId.index[f][z], globalId.index[altField(f)][y], *denseMat); // rows = time derivative of field ; cols = field vector multiplied 
				delete denseMat;
			}
		}

		if (opts_.fluxType == FluxType::Upwind) {


			#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
			#endif

			for (auto f : { E, H }) {
				auto intBdrFluxZeroNormal = buildByMult(*MInv[f], *buildZeroNormalIBFIOperator(f, model_, fes_, opts_), fes_);
				denseMat = intBdrFluxZeroNormal->SpMat().ToDenseMatrix();
				*denseMat *= -1.0;
				for (auto d : { X, Y, Z }) {
					globalOperator_.AddSubMatrix(globalId.index[f][d], globalId.index[f][d], *denseMat);
				}
				delete denseMat;
			}

			#ifdef SHOW_TIMER_INFORMATION
					std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
						(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
					std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
			#endif

			for (auto f : { E, H }) {
				for (auto d{ X }; d <= Z; d++) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						auto intBdrFluxTwoNormal = buildByMult(*MInv[f], *buildTwoNormalIBFIOperator(f, { d, d2 }, model_, fes_, opts_), fes_);
						auto denseMat = intBdrFluxTwoNormal->SpMat().ToDenseMatrix();
						globalOperator_.AddSubMatrix(globalId.index[f][d], globalId.index[f][d2], *denseMat); // rows = time derivative of field ; cols = field vector multiplied 
						delete denseMat;
					}
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto x{ X }; x <= Z; x++) { 
			auto y = (x + 1) % 3;
			auto z = (x + 2) % 3;
			auto dir = buildByMult(*MInv[f], *buildDerivativeOperator(x, fes_), fes_); // Example: field E, directional x, dt Ez, applied on Hy, +1.0.
			denseMat = dir->SpMat().ToDenseMatrix();
			*denseMat *= 1.0 - double(f) * 2.0; //cyclic permutations are positive for E, negative for H. Cyclic implies order: dt field polariz. (z) -> direction (x) -> applied field (y).
			globalOperator_.AddSubMatrix(globalId.index[f][z], globalId.index[altField(f)][y], *denseMat); // rows = time derivative of field ; cols = field vector multiplied 
			delete denseMat;
		}
		for (auto x{ X }; x <= Z; x++) {                                               // Example: field E, directional x, dt Ey, applied on Hz, +1.0.
			auto y = (x + 1) % 3;
			auto z = (x + 2) % 3;
			auto dir = buildByMult(*MInv[f], *buildDerivativeOperator(x, fes_), fes_);
			denseMat = dir->SpMat().ToDenseMatrix();
			*denseMat *= - 1.0 + double(f) * 2.0; //anticyclic permutations are negative for E, positive for H.
			globalOperator_.AddSubMatrix(globalId.index[f][y], globalId.index[altField(f)][z], *denseMat);
			delete denseMat;
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto x{ X }; x <= Z; x++) {
			auto y = (x + 1) % 3;
			auto z = (x + 2) % 3;
			auto fluxOneNormal = buildByMult(*MInv[f], *buildOneNormalOperator(altField(f), { x }, model_, fes_, opts_), fes_);
			denseMat = fluxOneNormal->SpMat().ToDenseMatrix();
			*denseMat *= 1.0 - double(f) * 2.0; //anticyclic permutations are positive for E, negative for H.
			globalOperator_.AddSubMatrix(globalId.index[f][y], globalId.index[altField(f)][z], *denseMat);
			delete denseMat;
		}
		for (auto x{ X }; x <= Z; x++) {
			auto y = (x + 1) % 3;
			auto z = (x + 2) % 3;
			auto fluxOneNormal = buildByMult(*MInv[f], *buildOneNormalOperator(altField(f), { x }, model_, fes_, opts_), fes_);
			denseMat = fluxOneNormal->SpMat().ToDenseMatrix();
			*denseMat *= - 1.0 + double(f) * 2.0; //cyclic permutations are negative for E, positive for H.
			globalOperator_.AddSubMatrix(globalId.index[f][z], globalId.index[altField(f)][y], *denseMat); // rows = time derivative of field ; cols = field vector multiplied 
			delete denseMat;
		}
	}

	if (opts_.fluxType == FluxType::Upwind) {

	#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
	#endif

		for (auto f : { E, H }) {
			auto fluxZeroNormal = buildByMult(*MInv[f], *buildZeroNormalOperator(f, model_, fes_, opts_), fes_);
			denseMat = fluxZeroNormal->SpMat().ToDenseMatrix();
			*denseMat *= -1.0;
			for (auto d : { X, Y, Z }) {
				globalOperator_.AddSubMatrix(globalId.index[f][d], globalId.index[f][d], *denseMat);
			}
			delete denseMat;
		}

	#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
	#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					auto fluxTwoNormal = buildByMult(*MInv[f], *buildTwoNormalOperator(f, { d, d2 }, model_, fes_, opts_), fes_);
					auto denseMat = fluxTwoNormal->SpMat().ToDenseMatrix();
					globalOperator_.AddSubMatrix(globalId.index[f][d], globalId.index[f][d2], *denseMat); // rows = time derivative of field; cols = field vector multiplied
					delete denseMat;
				}
			}
		}
    }

	globalOperator_.Finalize();

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Operator assembly finished" << std::endl;
	std::cout << std::endl;
#endif

	//CND_ = buildConductivityCoefficients(model_, fes_);

}

void MaxwellEvolution::Mult(const Vector& in, Vector& out) const
{
	const auto& dim{ fes_.GetMesh()->Dimension() };

	std::array<Vector, 3> eOld, hOld;
	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eNew[d].SetSpace(&fes_);
		hNew[d].SetSpace(&fes_);
		eNew[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hNew[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
		eNew[d] = 0.0;
		hNew[d] = 0.0;
	}

	globalOperator_.AddMult(in, out);

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<Planewave*>(source.get())) {
			
			auto func { evalTimeVarFunction(GetTime(),srcmngr_) };

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

				if (opts_.fluxType == FluxType::Upwind) {
					MFNN_GTFSF_[H][H][X][x]->Mult(func[H][X], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MFNN_GTFSF_[H][H][Y][x]->Mult(func[H][Y], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MFNN_GTFSF_[H][H][Z][x]->Mult(func[H][Z], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MP_GTFSF_[H]           ->Mult(func[H][x], hTemp[x]);
					hMap.TransferAdd(hTemp[x], hNew[x]);

					MFNN_GTFSF_[E][E][X][x]->Mult(func[E][X], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MFNN_GTFSF_[E][E][Y][x]->Mult(func[E][Y], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MFNN_GTFSF_[E][E][Z][x]->Mult(func[E][Z], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MP_GTFSF_[E]           ->Mult(func[E][x], eTemp[x]);
					eMap.TransferAdd(eTemp[x], eNew[x]);
				}

			}
		}
	}
}

}

