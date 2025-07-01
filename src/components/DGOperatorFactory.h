#pragma once

#include "ProblemDescription.h"
#include "evolution/GlobalEvolution.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "evolution/MaxwellEvolutionMethods.h"

#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"

#include "Types.h"

#include <chrono>
#include <iostream>
#include <string>

namespace maxwell
{

	using namespace mfem;
	using namespace mfemExtension;

	static FluxBdrCoefficientsCentered bdrCentCoeff{
		{BdrCond::PEC, {2.0, 0.0}},
		{BdrCond::PMC, {0.0, 2.0}},
		{BdrCond::SMA, {1.0, 1.0}},
		{BdrCond::SurfaceCond, {1.0, 1.0}},
	};

	static FluxBdrCoefficientsUpwind bdrUpwindCoeff{
		{BdrCond::PEC, {2.0, 0.0}},
		{BdrCond::PMC, {0.0, 2.0}},
		{BdrCond::SMA, {1.0, 1.0}},
		{BdrCond::SurfaceCond, {1.0, 1.0}},
	};

	static FluxSrcCoefficientsCentered srcCentCoeff{
		{BdrCond::TotalFieldIn, {1.0, 1.0}},
	};

	static FluxSrcCoefficientsUpwind srcUpwindCoeff{
		{BdrCond::TotalFieldIn, {1.0, 1.0}},
	};

	static FieldType altField(const FieldType &f)
	{
		switch (f)
		{
		case FieldType::E:
			return FieldType::H;
		case FieldType::H:
			return FieldType::E;
		default:
			throw std::runtime_error("Incorrect FieldType in input.");
		}
	}

	
	struct FieldOffsets
	{

		FieldOffsets(const int localBlockSize, const int nbrBlockSize, const FieldType f, const Direction d, const bool isLocal)
		{
			rowStartOffset = (3 * f + d) * localBlockSize;
			rowEndOffset = (3 * f + d) * localBlockSize + localBlockSize;
			colStartOffset = (3 * f + d) * (localBlockSize + nbrBlockSize);
			if (!isLocal){
				colEndOffset = (3 * f + d) * (localBlockSize + nbrBlockSize) + (localBlockSize + nbrBlockSize);
			}
			else{
				colEndOffset = (3 * f + d) * (localBlockSize + nbrBlockSize) + localBlockSize;
			}
		}

		int rowStartOffset;
		int rowEndOffset;
		int colStartOffset;
		int colEndOffset;
	};

	struct GlobalIndices
	{
		GlobalIndices(const int localBlockSize, const int nbrBlockSize, bool isLocal = false)
		{
			for (auto f : {E, H})
			{
				for (auto d : {X, Y, Z})
				{
					offsets[f][d] = std::make_unique<FieldOffsets>(localBlockSize, nbrBlockSize, f, d, isLocal);
				}
			}
		}

		std::array<std::array<std::unique_ptr<FieldOffsets>, 3>, 2> offsets;
	};

	static void loadBlockInGlobalAtIndices(const SparseMatrix &blk, SparseMatrix &dst, const std::pair<FieldOffsets, FieldOffsets> &ids, const double fieldSign)
	{
		auto expectedRows = ids.first.rowEndOffset - ids.first.rowStartOffset;
		auto expectedCols = ids.first.colEndOffset - ids.first.colStartOffset;
		MFEM_ASSERT(blk.NumRows() == expectedRows, "Block Sparse NumRows does not match intended number of Rows.");
		MFEM_ASSERT(blk.NumCols() == expectedCols, "Block Sparse NumCols does not match intended number of Cols.");
		Array<int> cols;
		Vector vals;
		for (auto r = 0; r < expectedRows; r++)
		{
			blk.GetRow(r, cols, vals);
			for (auto c = 0; c < cols.Size(); c++)
			{
				dst.Add(ids.first.rowStartOffset + r, ids.second.colStartOffset + cols[c], vals[c] * fieldSign);
			}
		}
	}

	static std::map<BdrCond, std::vector<double>> bdrCoeffCheck(double alpha)
	{
		std::map<BdrCond, std::vector<double>> res;
		if (alpha == 0.0)
		{
			res = bdrCentCoeff;
		}
		else
		{
			res = bdrUpwindCoeff;
		}
		return res;
	}

	template <typename FES, typename BF>
	std::unique_ptr<BF> buildByMult(
		const SparseMatrix &op1,
		const SparseMatrix &op2,
		FES &fes)
	{
		SparseMatrix *matrix = mfem::Mult(op1, op2);

		std::unique_ptr<BF> res = std::make_unique<BF>(&fes);

		res->Assemble();
		res->Finalize();
		res->SpMat().Swap(*matrix);
		delete matrix;

		return res;
	}

	void loadBlockInGlobalAtIndices(const SparseMatrix &blk, SparseMatrix &dst, const std::pair<Array<int>, Array<int>> &ids, const double fieldSign = 1.0);

	template <typename FES>
	class DGOperatorFactory
	{
	public:
		DGOperatorFactory(ProblemDescription &pd, FES &fes);

		// Methods for speficic FieldType or Direction Operators //
		template <typename BF>
		std::unique_ptr<BF> buildInverseMassMatrixSubOperator(const FieldType &f);

		template <typename BF>
		std::unique_ptr<BF> buildDerivativeSubOperator(const Direction &d);
		template <typename BF>
		std::unique_ptr<BF> buildZeroNormalSubOperator(const FieldType &f);
		template <typename BF>
		std::unique_ptr<BF> buildOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);
		template <typename BF>
		std::unique_ptr<BF> buildTwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);

		template <typename BF>
		std::unique_ptr<BF> buildZeroNormalIBFISubOperator(const FieldType &f);
		template <typename BF>
		std::unique_ptr<BF> buildOneNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);
		template <typename BF>
		std::unique_ptr<BF> buildTwoNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);

		// Methods for complete Maxwell Operators //
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellInverseMassMatrixOperator();
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellTFSFInverseMassMatrixOperator();

		template <typename BF>
		std::array<std::array<std::unique_ptr<BF>, 3>, 2> buildMaxwellDirectionalOperator();
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellZeroNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> buildMaxwellOneNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> buildMaxwellTwoNormalOperator();

		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellIntBdrZeroNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> buildMaxwellIntBdrOneNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> buildMaxwellIntBdrTwoNormalOperator();

		// Methors for complete Global Operators //
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildGlobalInverseMassMatrixOperator();

		template <typename BF>
		void addGlobalZeroNormalIBFIOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalOneNormalIBFIOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalTwoNormalIBFIOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalDirectionalOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalZeroNormalOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalOneNormalOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalTwoNormalOperators(mfem::SparseMatrix* global);

		std::unique_ptr<mfem::SparseMatrix> buildTFSFGlobalOperator();
		std::unique_ptr<mfem::SparseMatrix> buildGlobalOperator();

	private:
		ProblemDescription pd_;
		FES fes_;
	};

	template <typename FES>
	DGOperatorFactory<FES>::DGOperatorFactory(ProblemDescription &pd, FES &fes) : pd_(pd),
																				  fes_(fes)
	{
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildInverseMassMatrixSubOperator(const FieldType &f)
	{
		Vector aux{pd_.model.buildEpsMuPiecewiseVector(f)};
		PWConstCoefficient PWCoeff(aux);

		auto res = std::make_unique<BF>(&fes_);
		res->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildDerivativeSubOperator(const Direction &d)
	{
		auto res = std::make_unique<BF>(&fes_);

		if (d >= fes_.GetMesh()->Dimension())
		{
			res->Assemble();
			res->Finalize();
			return res;
		}

		ConstantCoefficient coeff = (d <= fes_.GetMesh()->Dimension()) ? ConstantCoefficient(1.0) : ConstantCoefficient(0.0);
		res->AddDomainIntegrator(
			new DerivativeIntegrator(coeff, d));

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildZeroNormalSubOperator(const FieldType &f)
	{
		auto res = std::make_unique<BF>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size())
		{
			for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
			{
				res->AddInteriorFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha), kv.second);
			}
		}
		else
		{
			res->AddInteriorFaceIntegrator(
				new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA)
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(c[kv.first].at(f) * pd_.opts.alpha), kv.second);
			}
			else
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(1.0), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size())
		{
			for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
			{
				res->AddInteriorFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), kv.second);
			}
		}
		else
		{
			res->AddInteriorFaceIntegrator(
				new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA)
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
			}
			else
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildTwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size())
		{
			for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
			{
				res->AddInteriorFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha), kv.second);
			}
		}
		else
		{
			res->AddInteriorFaceIntegrator(
				new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA)
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
			}
			else
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, 1.0), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildZeroNormalIBFISubOperator(const FieldType &f)
	{
		auto res = std::make_unique<BF>(&fes_);

		for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
		{
			if (kv.first != BdrCond::TotalFieldIn)
			{
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first)
				{
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator({}, 1.0), kv.second);
					break;
				default:
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator({}, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
					break;
				}
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildOneNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);

		for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
		{
			if (kv.first != BdrCond::TotalFieldIn)
			{
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first)
				{
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, 1.0), kv.second);
					break;
				default:
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
					break;
				}
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildTwoNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);

		for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
		{
			if (kv.first != BdrCond::TotalFieldIn)
			{
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first)
				{
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, 1.0), kv.second);
					break;
				default:
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
					break;
				}
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildMaxwellInverseMassMatrixOperator()
	{
		std::array<std::unique_ptr<BF>, 2> res;
		for (auto f : {E, H})
		{
			res[f] = buildInverseMassMatrixSubOperator<BF>(f);
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::unique_ptr<BF>, 3>, 2> DGOperatorFactory<FES>::buildMaxwellDirectionalOperator()
	{
		std::array<std::array<std::unique_ptr<BF>, 3>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				res[f][d] = buildByMult<FES, BF>(MInv->SpMat(), buildDerivativeSubOperator<BF>(d)->SpMat(), fes_);
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildMaxwellZeroNormalOperator()
	{
		std::array<std::unique_ptr<BF>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			res[f] = buildByMult<FES, BF>(MInv->SpMat(), buildZeroNormalSubOperator<BF>(f)->SpMat(), fes_);
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellOneNormalOperator()
	{
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					res[f][f2][d] = buildByMult<FES, BF>(MInv->SpMat(), buildOneNormalSubOperator<BF>(f2, {d})->SpMat(), fes_);
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellTwoNormalOperator()
	{
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					for (auto d2{X}; d2 <= Z; d2++)
					{
						res[f][f2][d][d2] = buildByMult<FES, BF>(MInv->SpMat(), buildTwoNormalSubOperator<BF>(f2, {d, d2})->SpMat(), fes_);
					}
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildMaxwellIntBdrZeroNormalOperator()
	{
		std::array<std::unique_ptr<BF>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			res[f] = buildByMult<FES, BF>(MInv->SpMat(), buildZeroNormalIBFISubOperator<BF>(f)->SpMat(), fes_);
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellIntBdrOneNormalOperator()
	{
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					res[f][f2][d] = buildByMult<FES, BF>(MInv->SpMat(), buildOneNormalIBFISubOperator<BF>(f2, {d})->SpMat(), fes_);
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellIntBdrTwoNormalOperator()
	{
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					for (auto d2{X}; d2 <= Z; d2++)
					{
						res[f][f2][d][d2] = buildByMult<FES, BF>(MInv->SpMat(), buildTwoNormalIBFISubOperator<BF>(f2, {d, d2})->SpMat(), fes_);
					}
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildGlobalInverseMassMatrixOperator()
	{
		return this->template  buildMaxwellInverseMassMatrixOperator<BF>();
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalZeroNormalIBFIOperators(SparseMatrix* global)
	{
		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		GlobalIndices globalId(fes_.GetNDofs(), additional_dofs);
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator<BF>();
			auto op = buildByMult<FES,BF>(
				MInv[f]->SpMat(), buildZeroNormalIBFISubOperator<BF>(f)->SpMat(), fes_);
			for (auto d : { X, Y, Z }) {
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d].get()),
					-1.0
				);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalOneNormalIBFIOperators(SparseMatrix* global)
	{
		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		GlobalIndices globalId(fes_.GetNDofs(), additional_dofs);
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator<BF>();
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(MInv[f]->SpMat(), buildOneNormalIBFISubOperator<BF>(altField(f), { x })->SpMat(), fes_);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][y].get(), *globalId.offsets[altField(f)][z].get()),
					1.0 - double(f) * 2.0
				);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][z].get(), *globalId.offsets[altField(f)][y].get()),
					-1.0 + double(f) * 2.0
				);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalTwoNormalIBFIOperators(SparseMatrix* global)
	{
		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		GlobalIndices globalId(fes_.GetNDofs(), additional_dofs);
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator<BF>();
			for (auto d{ X }; d <= Z; d++) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					auto op = buildByMult<FES,BF>(MInv[f]->SpMat(), buildTwoNormalIBFISubOperator<BF>(f, { d, d2 })->SpMat(), fes_);
					loadBlockInGlobalAtIndices(
						op->SpMat(),
						*global,
						std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d2].get()), 
						1.0
					);
				}
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalDirectionalOperators(SparseMatrix* global)
	{
		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		GlobalIndices globalId(fes_.GetNDofs(), additional_dofs, true);
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator<BF>();
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,ParBilinearForm>(
					MInv[f]->SpMat(), buildDerivativeSubOperator<BF>(x)->SpMat(), fes_);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][z].get(), *globalId.offsets[altField(f)][y].get()),
					1.0 - double(f) * 2.0
				);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][y].get(), *globalId.offsets[altField(f)][z].get()),
					-1.0 + double(f) * 2.0
				);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalZeroNormalOperators(SparseMatrix* global)
	{
		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		GlobalIndices globalId(fes_.GetNDofs(), additional_dofs);
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator<BF>();
			auto op = buildByMult<FES,BF>(
				MInv[f]->SpMat(), buildZeroNormalSubOperator<BF>(f)->SpMat(), fes_);
			for (auto d : { X, Y, Z }) {
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d].get()),
					-1.0
				);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalOneNormalOperators(SparseMatrix* global)
	{

		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		GlobalIndices globalId(fes_.GetNDofs(), additional_dofs);
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator<BF>();
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(
					MInv[f]->SpMat(), buildOneNormalSubOperator<BF>(altField(f), { x })->SpMat(), fes_);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][y].get(), *globalId.offsets[altField(f)][z].get()),
					1.0 - double(f) * 2.0
				);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(*globalId.offsets[f][z].get(), *globalId.offsets[altField(f)][y].get()),
					-1.0 + double(f) * 2.0
				);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalTwoNormalOperators(SparseMatrix* global)
	{

		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		GlobalIndices globalId(fes_.GetNDofs(), additional_dofs);
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator<BF>();
			for (auto d{ X }; d <= Z; d++) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					auto op = buildByMult<FES,BF>(
						MInv[f]->SpMat(), buildTwoNormalSubOperator<BF>(f, {d, d2})->SpMat(), fes_);
					loadBlockInGlobalAtIndices(
						op->SpMat(),
						*global,
						std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d2].get()),
						1.0
					);
				}
			}
		}
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildTFSFGlobalOperator()
	{

		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		std::unique_ptr<SparseMatrix> res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * (fes_.GetNDofs() + additional_dofs));

		std::chrono::high_resolution_clock::time_point startTime;
		#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
				startTime = std::chrono::high_resolution_clock::now();
				std::cout << "Assembling TFSF Inverse Mass One-Normal Operators" << std::endl;
		}
		#endif

		this->template addGlobalOneNormalOperators<BilinearForm>(res.get());

		#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
					std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
						(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
					startTime = std::chrono::high_resolution_clock::now();
					std::cout << "Assembling TFSF Inverse Mass Zero-Normal Operators" << std::endl;
		}
		#endif	

		this->template addGlobalZeroNormalOperators<BilinearForm>(res.get());


		#ifdef SHOW_TIMER_INFORMATION	
		if (Mpi::WorldRank() == 0){
					std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
						(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
					startTime = std::chrono::high_resolution_clock::now();
					std::cout << "Assembling TFSF Inverse Mass Two-Normal Operators" << std::endl;
		}
		#endif	

		this->template	addGlobalTwoNormalOperators<BilinearForm>(res.get());

		#ifdef SHOW_TIMER_INFORMATION	
		if (Mpi::WorldRank() == 0){
					std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
						(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		}
		#endif

		res->Finalize();
		return res;
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildGlobalOperator()
	{

		auto additional_dofs = 0;
		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
        	additional_dofs = fes_.num_face_nbr_dofs;
    	}

		std::unique_ptr<SparseMatrix> res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * (fes_.GetNDofs() + additional_dofs));


		if (pd_.model.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

			std::chrono::high_resolution_clock::time_point startTime;
			#ifdef SHOW_TIMER_INFORMATION
			if (Mpi::WorldRank() == 0){
				startTime = std::chrono::high_resolution_clock::now() ;
			}
			#endif

			#ifdef SHOW_TIMER_INFORMATION
			if (Mpi::WorldRank() == 0){
						std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
			}
			#endif

				this->template	addGlobalOneNormalIBFIOperators<ParBilinearForm>(res.get());

			#ifdef SHOW_TIMER_INFORMATION
			if (Mpi::WorldRank() == 0){
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
				startTime = std::chrono::high_resolution_clock::now();
				std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
			}
			#endif

				this->template	addGlobalZeroNormalIBFIOperators<ParBilinearForm>(res.get());

			#ifdef SHOW_TIMER_INFORMATION
			if (Mpi::WorldRank() == 0){
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
				startTime = std::chrono::high_resolution_clock::now();
				std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
			}
			#endif

				this->template	addGlobalTwoNormalIBFIOperators<ParBilinearForm>(res.get());

			#ifdef SHOW_TIMER_INFORMATION
			if (Mpi::WorldRank() == 0){
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			}
			#endif

		}
		else{

			#ifdef SHOW_TIMER_INFORMATION		
					if (Mpi::WorldRank() == 0){
							std::cout << "No Interior Boundary Operators to Assemble." << std::endl;
					}
			#endif

		}

		std::chrono::high_resolution_clock::time_point startTime;
		#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
		}
		#endif

		this->template	addGlobalDirectionalOperators<ParBilinearForm>(res.get());

		#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
		}
		#endif

		this->template	addGlobalOneNormalOperators<ParBilinearForm>(res.get());

		#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
		}
		#endif

		this->template	addGlobalZeroNormalOperators<ParBilinearForm>(res.get());

		#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
		}
		#endif

		this->template	addGlobalTwoNormalOperators<ParBilinearForm>(res.get());

		#ifdef SHOW_TIMER_INFORMATION
		if (Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "All operators assembled. Finalising global operator matrix." << std::endl;
		}
		#endif

		res->Finalize();
		return res;
	}

}