#include "Evolution.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

Evolution::Evolution(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{

	for (auto f : { E, H }) {
		MP_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildZeroNormalOperator(f, model_, fes_, opts_), fes_);
		for (auto d{ X }; d <= Z; d++) {
			MS_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(d, fes_), fes_);
			for (auto d2{ X }; d2 <= Z; d2++) {
				for (auto f2 : { E, H }) {
					MFN_[f][f2][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildOneNormalOperator(f2, { d }, model_, fes_, opts_), fes_);
					MFNN_[f][f2][d][d2] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildTwoNormalOperator(f2, { d, d2 }, model_, fes_, opts_), fes_);
				}
			}
		}
	}

	for (auto f : { E, H }) {//TFSF - SrcConds
		MFF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyFixOperator(f, {}, model_, fes_, opts_), fes_);
		for (auto d : { X, Y, Z }) {
			MBF_[f][d] = buildByMult(
				*buildInverseMassMatrix(f, model_, fes_), *buildFluxFunctionOperator(f, { X }, model_, fes_, opts_), fes_);
		}
	}

	for (auto bdr_att = 0; bdr_att < model_.getConstMesh().GetNBE(); bdr_att++) {
		if (model_.getConstMesh().GetBdrAttribute(bdr_att) == 301) {
			srcmngr_.initTFSFPreReqs(model_.getConstMesh());
			auto fesTF{ srcmngr_.getTFSpace() };
			auto fesSF{ srcmngr_.getSFSpace() };
			
			for (auto f : { E, H }) {
				MP_TF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, *fesTF), *buildZeroNormalOperator(f, model_, *fesTF, opts_), *fesTF);
				for (auto d{ X }; d <= Z; d++) {
					MS_TF_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, *fesTF), *buildDerivativeOperator(d, *fesTF), *fesTF);
					for (auto d2{ X }; d2 <= Z; d2++) {
						for (auto f2 : { E, H }) {
							MFN_TF_[f][f2][d] = buildByMult(*buildInverseMassMatrix(f, model_, *fesTF), *buildOneNormalOperator(f2, { d }, model_, *fesTF, opts_), *fesTF);
							MFNN_TF_[f][f2][d][d2] = buildByMult(*buildInverseMassMatrix(f, model_, *fesTF), *buildTwoNormalOperator(f2, { d, d2 }, model_, *fesTF, opts_), *fesTF);
						}
					}
				}
			}

			for (auto f : { E, H }) {
				MP_SF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, *fesSF), *buildZeroNormalOperator(f, model_, *fesSF, opts_), *fesSF);
				for (auto d{ X }; d <= Z; d++) {
					MS_SF_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, *fesSF), *buildDerivativeOperator(d, *fesSF), *fesSF);
					for (auto d2{ X }; d2 <= Z; d2++) {
						for (auto f2 : { E, H }) {
							MFN_SF_[f][f2][d] = buildByMult(*buildInverseMassMatrix(f, model_, *fesSF), *buildOneNormalOperator(f2, { d }, model_, *fesSF, opts_), *fesSF);
							MFNN_SF_[f][f2][d][d2] = buildByMult(*buildInverseMassMatrix(f, model_, *fesSF), *buildTwoNormalOperator(f2, { d, d2 }, model_, *fesSF, opts_), *fesSF);
						}
					}
				}
			}

			for (auto f : { E, H }) {
				MTF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, *fesTF), *buildTFSFOperator(f, *fesTF,  1.0), *fesTF);
				MSF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, *fesSF), *buildTFSFOperator(f, *fesSF, -1.0), *fesSF);
			}

			break;
		}
	}

	if (model_.getInteriorBoundaryToMarker().size() != 0) {
		for (auto f : { E, H }) { //IntBdrConds
			MPB_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildZeroNormalIBFIOperator(f, model_, fes_, opts_), fes_);
			for (auto d{ X }; d <= Z; d++) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					for (auto f2 : { E, H }) {
						MFNB_[f][f2][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildOneNormalIBFIOperator(f2, { d }, model_, fes_, opts_), fes_);
						MFNNB_[f][f2][d][d2] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildTwoNormalIBFIOperator(f2, { d, d2 }, model_, fes_, opts_), fes_);
					}
				}
			}
		}
	}

 }

void Evolution::Mult(const Vector& in, Vector& out) const
{
	const auto& dim{ fes_.GetMesh()->Dimension() };

	std::array<Vector, 3> eOld, hOld;
	std::array<GridFunction, 3> eNew, hNew;
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
		MS_[H][y]		 ->AddMult(eOld[z], hNew[x],-1.0);
		MS_[H][z]		 ->AddMult(eOld[y], hNew[x]);
		MS_[E][y]		 ->AddMult(hOld[z], eNew[x]);
		MS_[E][z]		 ->AddMult(hOld[y], eNew[x],-1.0);
		
		MFN_[H][E][y]  	 ->AddMult(eOld[z], hNew[x], 1.0);
		MFN_[H][E][z]	 ->AddMult(eOld[y], hNew[x],-1.0);
		MFN_[E][H][y]  	 ->AddMult(hOld[z], eNew[x],-1.0);
		MFN_[E][H][z]	 ->AddMult(hOld[y], eNew[x], 1.0);

		if (opts_.fluxType == FluxType::Upwind) {

			MFNN_[H][H][X][x]->AddMult(hOld[X], hNew[x], 1.0);
			MFNN_[H][H][Y][x]->AddMult(hOld[Y], hNew[x], 1.0);
			MFNN_[H][H][Z][x]->AddMult(hOld[Z], hNew[x], 1.0);
			MP_[H]			 ->AddMult(hOld[x], hNew[x],-1.0);
			
			MFNN_[E][E][Y][x]->AddMult(eOld[Y], eNew[x], 1.0);
			MFNN_[E][E][X][x]->AddMult(eOld[X], eNew[x], 1.0);
			MFNN_[E][E][Z][x]->AddMult(eOld[Z], eNew[x], 1.0);
			MP_[E]			 ->AddMult(eOld[x], eNew[x],-1.0);
		}

		if (model_.getInteriorBoundaryToMarker().size() != 0) {
			
			MFNB_[H][E][y]->AddMult(eOld[z], hNew[x]);
			MFNB_[H][E][z]->AddMult(eOld[y], hNew[x], -1.0);
			MFNB_[E][H][y]->AddMult(hOld[z], eNew[x], -1.0);
			MFNB_[E][H][z]->AddMult(hOld[y], eNew[x]);

			if (opts_.fluxType == FluxType::Upwind) {

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

	}

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<Planewave*>(source.get())) {
			
			auto time{ GetTime() };
			
			auto func_tf{ srcmngr_.evalTimeVarField(time, true) };
			auto func_sf{ srcmngr_.evalTimeVarField(time, false) };
			
			std::array<GridFunction, 3> eTempTF, hTempTF, eTempSF, hTempSF;
			
			for (int d = X; d <= Z; d++) {
				eTempTF[d].SetSpace(srcmngr_.getTFSpace());
				hTempTF[d].SetSpace(srcmngr_.getTFSpace());
				eTempSF[d].SetSpace(srcmngr_.getSFSpace());
				hTempSF[d].SetSpace(srcmngr_.getSFSpace());
				eTempTF[d] = 0.0;
				hTempTF[d] = 0.0;
				eTempSF[d] = 0.0;
				hTempSF[d] = 0.0;
			}

			for (int x = X; x <= Z; x++) {

				MTF_[E]->Mult(func_tf[E][x], eTempTF[x]);
				MTF_[H]->Mult(func_tf[H][x], hTempTF[x]);

				MSF_[E]->Mult(func_sf[E][x], eTempSF[x]);
				MSF_[H]->Mult(func_sf[H][x], hTempSF[x]);

				MaxwellTransferMap eMapTF(eTempTF[x], eNew[x]);
				eMapTF.TransferAdd		 (eTempTF[x], eNew[x]);

				MaxwellTransferMap hMapTF(hTempTF[x], hNew[x]);
				hMapTF.TransferAdd		 (hTempTF[x], hNew[x]);

				MaxwellTransferMap eMapSF(eTempSF[x], eNew[x]);
				eMapSF.TransferAdd		 (eTempSF[x], eNew[x]);

				MaxwellTransferMap hMapSF(hTempSF[x], hNew[x]);
				hMapSF.TransferAdd		 (hTempSF[x], hNew[x]);

			}
		}
	}

}

}

