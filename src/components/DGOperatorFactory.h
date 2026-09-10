#pragma once

#include "ProblemDescription.h"
#include "SCPMLLayout.h"

#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <filesystem>
#include <type_traits>
#include <vector>

namespace maxwell
{

	using namespace mfem;
	using namespace mfemExtension;

	inline const FluxBdrCoefficientsCentered bdrCentCoeff{
		{BdrCond::PEC, {2.0, 0.0}},
		{BdrCond::PMC, {0.0, 2.0}},
		{BdrCond::SMA, {1.0, 1.0}},
		{BdrCond::SurfaceCond, {1.0, 1.0}},
		{BdrCond::SGBC, {1.0, 1.0}}
	};

	inline const FluxBdrCoefficientsUpwind bdrUpwindCoeff{
		{BdrCond::PEC, {2.0, 0.0}},
		{BdrCond::PMC, {0.0, 2.0}},
		{BdrCond::SMA, {1.0, 1.0}},
		{BdrCond::SurfaceCond, {1.0, 1.0}},
		{BdrCond::SGBC, {1.0, 1.0}}
	};

	inline const FluxSrcCoefficientsCentered srcCentCoeff{
		{BdrCond::TotalFieldIn, {1.0, 1.0}},
		{BdrCond::SGBC, {1.0, 1.0}},
	};

	inline const FluxSrcCoefficientsUpwind srcUpwindCoeff{
		{BdrCond::TotalFieldIn, {1.0, 1.0}},
		{BdrCond::SGBC, {1.0, 1.0}},
	};

	inline FieldType altField(const FieldType &f)
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

	inline void loadBlockInGlobalAtIndices(const SparseMatrix &blk, SparseMatrix &dst, const std::pair<FieldOffsets, FieldOffsets> &ids, const double fieldSign)
	{
		auto expectedRows = ids.first.rowEndOffset - ids.first.rowStartOffset;
		auto expectedCols = ids.first.colEndOffset - ids.first.colStartOffset;
		MFEM_ASSERT(blk.NumRows() == expectedRows, "Block Sparse NumRows does not match intended number of Rows.");
		MFEM_ASSERT(blk.NumCols() >= expectedCols, "Block Sparse NumCols is smaller than intended number of Cols.");
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

	inline std::map<BdrCond, std::vector<double>> bdrCoeffCheck(double alpha)
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

	/// A block placement for CSR-direct assembly (S1 optimization).
	/// Stores a finalized sub-operator and its position within the global matrix.
	struct CSRBlockPlacement {
		std::unique_ptr<SparseMatrix> block;  ///< Finalized CSR sub-operator.
		int rowOffset;                        ///< Global row offset for this block.
		int colOffset;                        ///< Global column offset for this block.
		double sign;                          ///< Scaling factor (±1).
	};

	/// Collect a sub-operator block placement instead of scattering into a LIL matrix.
	/// The block's CSR data is copied (for operators placed at multiple offsets).
	inline void collectBlockPlacement(
		const SparseMatrix& blk,
		std::vector<CSRBlockPlacement>& blocks,
		const std::pair<FieldOffsets, FieldOffsets>& ids,
		double fieldSign)
	{
		blocks.push_back(CSRBlockPlacement{
			std::make_unique<SparseMatrix>(blk),
			ids.first.rowStartOffset,
			ids.second.colStartOffset,
			fieldSign
		});
	}

	inline void collectBlockPlacement(
		const SparseMatrix& blk,
		std::vector<CSRBlockPlacement>& blocks,
		int rowOffset,
		int colOffset,
		double fieldSign)
	{
		blocks.push_back(CSRBlockPlacement{
			std::make_unique<SparseMatrix>(blk),
			rowOffset,
			colOffset,
			fieldSign
		});
	}

	/// Merge collected CSR block placements into a single finalized CSR SparseMatrix.
	/// Uses a two-pass marker technique (like MFEM's Add) to avoid LIL overhead.
	inline std::unique_ptr<SparseMatrix> mergeBlocksToCSR(
		std::vector<CSRBlockPlacement>& blocks,
		int globalRows, int globalCols)
	{
		// Pass 1: Count unique column entries per global row.
		std::vector<int> marker(globalCols, -1);
		std::vector<int64_t> row_ptr(static_cast<size_t>(globalRows) + 1, 0);

		for (int row = 0; row < globalRows; ++row) {
			int nnz = 0;
			for (auto& bp : blocks) {
				const int localRow = row - bp.rowOffset;
				if (localRow < 0 || localRow >= bp.block->Height()) continue;
				const int rowNNZ = bp.block->RowSize(localRow);
				const int* cols = bp.block->GetRowColumns(localRow);
				for (int k = 0; k < rowNNZ; ++k) {
					int globalCol = cols[k] + bp.colOffset;
					if (marker[globalCol] != row) {
						marker[globalCol] = row;
						++nnz;
					}
				}
			}
			row_ptr[row + 1] = row_ptr[row] + nnz;
		}

		const int64_t totalNNZ_64 = row_ptr[globalRows];
		if (totalNNZ_64 <= 0 ||
		    totalNNZ_64 > static_cast<int64_t>(std::numeric_limits<int>::max())) {
			int block_nnz = 0;
			for (const auto& bp : blocks) {
				if (bp.block) {
					block_nnz += bp.block->NumNonZeroElems();
				}
			}
			throw std::runtime_error(
				"mergeBlocksToCSR: total NNZ (" + std::to_string(totalNNZ_64) +
				") is out of range for SparseMatrix indexing."
				" rank=" + std::to_string(Mpi::WorldRank()) +
				" nblocks=" + std::to_string(blocks.size()) +
				" sum_block_nnz=" + std::to_string(block_nnz) +
				" rows=" + std::to_string(globalRows) +
				" cols=" + std::to_string(globalCols));
		}
		const int totalNNZ = static_cast<int>(totalNNZ_64);

		std::vector<int> C_i(static_cast<size_t>(globalRows) + 1);
		for (int i = 0; i <= globalRows; ++i) {
			C_i[i] = static_cast<int>(row_ptr[i]);
		}

		std::vector<int> C_j(static_cast<size_t>(totalNNZ));
		std::vector<real_t> C_data(static_cast<size_t>(totalNNZ));

		// Pass 2: Fill J and A arrays, merging duplicate column entries.
		std::fill(marker.begin(), marker.end(), -1);
		int pos = 0;
		for (int row = 0; row < globalRows; ++row) {
			for (auto& bp : blocks) {
				const int localRow = row - bp.rowOffset;
				if (localRow < 0 || localRow >= bp.block->Height()) continue;
				const int rowNNZ = bp.block->RowSize(localRow);
				const int* cols = bp.block->GetRowColumns(localRow);
				const real_t* vals = bp.block->GetRowEntries(localRow);
				for (int k = 0; k < rowNNZ; ++k) {
					int globalCol = cols[k] + bp.colOffset;
					if (marker[globalCol] < C_i[row]) {
						// New entry for this row.
						C_j[pos] = globalCol;
						C_data[pos] = vals[k] * bp.sign;
						marker[globalCol] = pos;
						++pos;
					} else {
						// Duplicate column — accumulate.
						C_data[marker[globalCol]] += vals[k] * bp.sign;
					}
				}
			}
		}

		int* C_i_ptr = mfem::Memory<int>(globalRows + 1);
		int* C_j_ptr = mfem::Memory<int>(totalNNZ);
		real_t* C_data_ptr = mfem::Memory<real_t>(totalNNZ);
		std::memcpy(C_i_ptr, C_i.data(), static_cast<size_t>(globalRows + 1) * sizeof(int));
		std::memcpy(C_j_ptr, C_j.data(), static_cast<size_t>(totalNNZ) * sizeof(int));
		std::memcpy(C_data_ptr, C_data.data(), static_cast<size_t>(totalNNZ) * sizeof(real_t));

		return std::make_unique<SparseMatrix>(C_i_ptr, C_j_ptr, C_data_ptr, globalRows, globalCols);
	}

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

		template <typename BF>
		std::unique_ptr<BF> buildSourceFaceIBFIZeroNormalSubOperator(const FieldType &f, mfem::Array<int>& marker);
		template <typename BF>
		std::unique_ptr<BF> buildSourceFaceIBFIOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker);
		template <typename BF>
		std::unique_ptr<BF> buildSourceFaceIBFITwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker);

		template <typename BF>
		std::unique_ptr<BF> buildBoundarySourceFaceIBFIZeroNormalSubOperator(const FieldType &f, mfem::Array<int>& marker);
		template <typename BF>
		std::unique_ptr<BF> buildBoundarySourceFaceIBFIOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker);
		template <typename BF>
		std::unique_ptr<BF> buildBoundarySourceFaceIBFITwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker);

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
		std::unique_ptr<BF> buildSigmaMassOperator();

		template <typename BF>
		void addGlobalZeroNormalIBFIOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalOneNormalIBFIOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalTwoNormalIBFIOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalSourceFaceIBFIZeroNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker);
		template <typename BF>
		void addGlobalSourceFaceIBFIOneNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker);
		template <typename BF>
		void addGlobalSourceFaceIBFITwoNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker);
		template <typename BF>
		void addGlobalDirectionalOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalZeroNormalOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalOneNormalOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalTwoNormalOperators(mfem::SparseMatrix* global);
		template <typename BF>
		void addGlobalConductiveOperator(mfem::SparseMatrix* global);

		// Overloads accepting pre-computed M^{-1} to avoid redundant rebuilds.
		template <typename BF>
		void addGlobalZeroNormalIBFIOperators(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalOneNormalIBFIOperators(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalTwoNormalIBFIOperators(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalDirectionalOperators(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalZeroNormalOperators(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalOneNormalOperators(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalTwoNormalOperators(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalConductiveOperator(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalSourceFaceIBFIZeroNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalSourceFaceIBFIOneNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalSourceFaceIBFITwoNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv); 
		template <typename BF>
		void addGlobalBoundarySourceFaceIBFIZeroNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalBoundarySourceFaceIBFIOneNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void addGlobalBoundarySourceFaceIBFITwoNormalOperators(mfem::SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv);

		// S1: Overloads that collect block placements for CSR-direct assembly.
		template <typename BF>
		void collectGlobalZeroNormalIBFIOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void collectGlobalOneNormalIBFIOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void collectGlobalTwoNormalIBFIOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void collectGlobalDirectionalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void collectGlobalZeroNormalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void collectGlobalOneNormalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void collectGlobalTwoNormalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);
		template <typename BF>
		void collectGlobalConductiveOperator(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv);

		std::unique_ptr<mfem::SparseMatrix> buildTFSFGlobalOperator();
		std::unique_ptr<mfem::SparseMatrix> buildSGBCGlobalOperator();
		std::unique_ptr<mfem::SparseMatrix> buildSourceFaceOperator(BdrCond filter);
		std::unique_ptr<mfem::SparseMatrix> buildSourceFaceOperator(mfem::Array<int>& marker);
		std::unique_ptr<mfem::SparseMatrix> buildGlobalOperator();
		/// Bagci/Chen SC-PML ADE + optional curl a-rescale (κ>1).
		/// curl_delta[u] is ndofs×ndofs: out_Fu += Delta_u * out_Fu for F in {E,H}.
		/// Entries are null when all regions have kappa_max == 1.
		void buildSCPMLOperators(
			const SCPMLLayout& layout,
			std::unique_ptr<mfem::SparseMatrix>& ade_operator,
			std::array<std::unique_ptr<mfem::SparseMatrix>, 3>& curl_delta);

	private:
		ProblemDescription pd_;
		FES fes_;

		template <typename BF>
		std::unique_ptr<BF> buildMarkedMassOperator(
			mfem::Coefficient& coeff, mfem::Array<int>& attr_marker);

		template <typename BF>
		std::unique_ptr<BF> buildMarkedInverseMassOperator(
			mfem::Coefficient& coeff, mfem::Array<int>& attr_marker);

		mfem::Array<int> buildInteriorIgnoreMarker() const
		{
			mfem::Array<int> marker;
			if (!pd_.model.getInteriorBoundaryToMarker().empty()) {
				int marker_size = pd_.model.getInteriorBoundaryToMarker().begin()->second.Size();
				marker.SetSize(marker_size);
				marker = 0;
				for (const auto &kv : pd_.model.getInteriorBoundaryToMarker()) {
					if (kv.first != BdrCond::TotalFieldIn) {
						for (int i = 0; i < kv.second.Size(); i++) {
							if (kv.second[i] == 1) marker[i] = 1;
						}
					}
				}
			}
			return marker;
		}

		int getAdditionalDofs() const
		{
			if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
				return fes_.num_face_nbr_dofs;
			}
			return 0;
		}

		int meshDimension() const
		{
			return fes_.GetMesh()->Dimension();
		}
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
		auto* integ = new DerivativeIntegrator(coeff, d);

		// For curved (high-order geometry) meshes, MFEM's default DerivativeIntegrator
		// quadrature rule (order 2p-1 for Pk) does not account for the non-constant
		// Jacobian. The adjugate matrix adj(J) introduces extra polynomial degree
		// (dim-1)*(meshOrder-1). Under-integration breaks the discrete summation-by-parts
		// (SBP) property, causing a slow-growing energy instability on curved elements.
		auto* nodalFES = fes_.GetMesh()->GetNodalFESpace();
		if (nodalFES && fes_.GetMesh()->GetNE() > 0) {
			int meshOrder = nodalFES->GetMaxElementOrder();
			if (meshOrder > 1) {
				int p = fes_.FEColl()->GetOrder();
				int dim = fes_.GetMesh()->Dimension();
				int adjDeg = (dim - 1) * (meshOrder - 1);
				int totalOrder = 2 * p - 1 + adjDeg;
				auto geomType = fes_.GetMesh()->GetElementGeometry(0);
				integ->SetIntRule(&IntRules.Get(geomType, totalOrder));
			}
		}

		res->AddDomainIntegrator(integ);

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildZeroNormalSubOperator(const FieldType &f)
	{
		auto res = std::make_unique<BF>(&fes_);

		auto ignore_marker = buildInteriorIgnoreMarker();

		if (ignore_marker.Size() > 0) {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha), ignore_marker);
		} else {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{
			if (kv.first == BdrCond::SGBC)
			{
				continue;
			}
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

		auto ignore_marker = buildInteriorIgnoreMarker();

		if (ignore_marker.Size() > 0) {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), ignore_marker);
		} else {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{
			if (kv.first == BdrCond::SGBC)
			{
				continue;
			}
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

		auto ignore_marker = buildInteriorIgnoreMarker();

		if (ignore_marker.Size() > 0) {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha), ignore_marker);
		} else {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{
			if (kv.first == BdrCond::SGBC)
			{
				continue;
			}
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
            if (kv.first != BdrCond::TotalFieldIn && kv.first != BdrCond::SGBC)
            {
                auto c = bdrCoeffCheck(pd_.opts.alpha);
                switch (kv.first) {
                    case BdrCond::SMA:
                    case BdrCond::PEC:
                    case BdrCond::PMC:
                    {
                        double coeff = (kv.first == BdrCond::SMA) ? 1.0 : (c[kv.first].at(f) * pd_.opts.alpha);
                        res->AddInternalBoundaryFaceIntegrator(
                            new mfemExtension::MaxwellDGInteriorJumpIntegrator({}, coeff), kv.second);
                        break;
                    }
                    default:
                        res->AddInternalBoundaryFaceIntegrator(
                            new mfemExtension::MaxwellDGZeroNormalJumpIntegrator(c[kv.first].at(f) * pd_.opts.alpha), kv.second);
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
            if (kv.first != BdrCond::TotalFieldIn && kv.first != BdrCond::SGBC)
            {
                auto c = bdrCoeffCheck(pd_.opts.alpha);
                switch (kv.first) {
                    case BdrCond::SMA:
                    case BdrCond::PEC:
                    case BdrCond::PMC:
                    {
                        double coeff = (kv.first == BdrCond::SMA) ? 1.0 : c[kv.first].at(f);
                        res->AddInternalBoundaryFaceIntegrator(
                            new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, coeff), kv.second);
                        break;
                    }
                    default:
                        res->AddInternalBoundaryFaceIntegrator(
                            new mfemExtension::MaxwellDGOneNormalJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
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
            if (kv.first != BdrCond::TotalFieldIn && kv.first != BdrCond::SGBC)
            {
                auto c = bdrCoeffCheck(pd_.opts.alpha);
                switch (kv.first) {
                    case BdrCond::SMA:
                    case BdrCond::PEC:
                    case BdrCond::PMC:
                    {
                        double coeff = (kv.first == BdrCond::SMA) ? 1.0 : (c[kv.first].at(f) * pd_.opts.alpha);
                        res->AddInternalBoundaryFaceIntegrator(
                            new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, coeff), kv.second);
                        break;
                    }
                    default:
                        res->AddInternalBoundaryFaceIntegrator(
                            new mfemExtension::MaxwellDGTwoNormalJumpIntegrator(dirTerms, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
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
    std::unique_ptr<BF> DGOperatorFactory<FES>::buildSourceFaceIBFIZeroNormalSubOperator(const FieldType &f, mfem::Array<int>& marker)
    {
        auto res = std::make_unique<BF>(&fes_);
        res->AddInternalBoundaryFaceIntegrator(
            new mfemExtension::MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha), marker);
        res->Assemble();
        res->Finalize();
        return res;
    }

    template <typename FES>
    template <typename BF>
    std::unique_ptr<BF> DGOperatorFactory<FES>::buildSourceFaceIBFIOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker)
    {
        auto res = std::make_unique<BF>(&fes_);
        res->AddInternalBoundaryFaceIntegrator(
            new mfemExtension::MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), marker);
        res->Assemble();
        res->Finalize();
        return res;
    }

    template <typename FES>
    template <typename BF>
    std::unique_ptr<BF> DGOperatorFactory<FES>::buildSourceFaceIBFITwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker)
    {
        auto res = std::make_unique<BF>(&fes_);
        res->AddInternalBoundaryFaceIntegrator(
            new mfemExtension::MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha), marker);
        res->Assemble();
        res->Finalize();
        return res;
    }

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildBoundarySourceFaceIBFIZeroNormalSubOperator(const FieldType &f, mfem::Array<int>& marker)
	{
		auto res = std::make_unique<BF>(&fes_);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha), marker);
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildBoundarySourceFaceIBFIOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker)
	{
		auto res = std::make_unique<BF>(&fes_);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), marker);
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildBoundarySourceFaceIBFITwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms, mfem::Array<int>& marker)
	{
		auto res = std::make_unique<BF>(&fes_);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha), marker);
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
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildSigmaMassOperator()
	{
		Vector sigma = pd_.model.buildSigmaPiecewiseVector(); 
		PWConstCoefficient SigCoeff(sigma);

		auto bf = std::make_unique<BF>(&fes_);
		bf->AddDomainIntegrator(new MassIntegrator(SigCoeff));
		bf->Assemble();
		bf->Finalize();
		return bf;
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalZeroNormalIBFIOperators(SparseMatrix* global)
	{
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalZeroNormalIBFIOperators<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalZeroNormalIBFIOperators(SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
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
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalOneNormalIBFIOperators<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalOneNormalIBFIOperators(SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue; // S2: normal component x is zero in lower dimensions
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
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalTwoNormalIBFIOperators<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalTwoNormalIBFIOperators(SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				if (d >= dim) continue; // S2: normal component d is zero in lower dimensions
				for (auto d2{ X }; d2 <= Z; d2++) {
					if (d2 >= dim) continue; // S2: normal component d2 is zero in lower dimensions
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
	void DGOperatorFactory<FES>::addGlobalSourceFaceIBFIZeroNormalOperators(SparseMatrix* global, mfem::Array<int>& marker)
	{
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalSourceFaceIBFIZeroNormalOperators<BF>(global, marker, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalSourceFaceIBFIZeroNormalOperators(SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			auto op = buildByMult<FES,BF>(
				MInv[f]->SpMat(), buildSourceFaceIBFIZeroNormalSubOperator<BF>(f, marker)->SpMat(), fes_);
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
	void DGOperatorFactory<FES>::addGlobalSourceFaceIBFIOneNormalOperators(SparseMatrix* global, mfem::Array<int>& marker)
	{
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalSourceFaceIBFIOneNormalOperators<BF>(global, marker, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalSourceFaceIBFIOneNormalOperators(SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue;
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(MInv[f]->SpMat(), buildSourceFaceIBFIOneNormalSubOperator<BF>(altField(f), { x }, marker)->SpMat(), fes_);
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
	void DGOperatorFactory<FES>::addGlobalSourceFaceIBFITwoNormalOperators(SparseMatrix* global, mfem::Array<int>& marker)
	{
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalSourceFaceIBFITwoNormalOperators<BF>(global, marker, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalSourceFaceIBFITwoNormalOperators(SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				if (d >= dim) continue;
				for (auto d2{ X }; d2 <= Z; d2++) {
					if (d2 >= dim) continue;
					auto op = buildByMult<FES,BF>(MInv[f]->SpMat(), buildSourceFaceIBFITwoNormalSubOperator<BF>(f, { d, d2 }, marker)->SpMat(), fes_);
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
	void DGOperatorFactory<FES>::addGlobalBoundarySourceFaceIBFIZeroNormalOperators(SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			auto op = buildByMult<FES,BF>(
				MInv[f]->SpMat(), buildBoundarySourceFaceIBFIZeroNormalSubOperator<BF>(f, marker)->SpMat(), fes_);
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
	void DGOperatorFactory<FES>::addGlobalBoundarySourceFaceIBFIOneNormalOperators(SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue;
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(
					MInv[f]->SpMat(), buildBoundarySourceFaceIBFIOneNormalSubOperator<BF>(altField(f), { x }, marker)->SpMat(), fes_);
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
	void DGOperatorFactory<FES>::addGlobalBoundarySourceFaceIBFITwoNormalOperators(SparseMatrix* global, mfem::Array<int>& marker, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				if (d >= dim) continue;
				for (auto d2{ X }; d2 <= Z; d2++) {
					if (d2 >= dim) continue;
					auto op = buildByMult<FES,BF>(
						MInv[f]->SpMat(), buildBoundarySourceFaceIBFITwoNormalSubOperator<BF>(f, { d, d2 }, marker)->SpMat(), fes_);
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
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalDirectionalOperators<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalDirectionalOperators(SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue; // S2: derivative in direction x is zero beyond mesh dimension
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(
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
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalZeroNormalOperators<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalZeroNormalOperators(SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
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
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalOneNormalOperators<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalOneNormalOperators(SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue; // S2: normal component x is zero in lower dimensions
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
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalTwoNormalOperators<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalTwoNormalOperators(SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				if (d >= dim) continue; // S2: normal component d is zero in lower dimensions
				for (auto d2{ X }; d2 <= Z; d2++) {
					if (d2 >= dim) continue; // S2: normal component d2 is zero in lower dimensions
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
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalConductiveOperator(mfem::SparseMatrix* global)
	{
		auto MInv = buildMaxwellInverseMassMatrixOperator<BF>();
		addGlobalConductiveOperator<BF>(global, MInv);
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::addGlobalConductiveOperator(mfem::SparseMatrix* global, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		auto MSig  = buildSigmaMassOperator<BF>();
		auto ASigE = buildByMult<FES, BF>(MInv[E]->SpMat(), MSig->SpMat(), fes_);

		GlobalIndices gid(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto d : { X, Y, Z }) {
			loadBlockInGlobalAtIndices(
				ASigE->SpMat(),
				*global,
				std::make_pair(*gid.offsets[E][d].get(), *gid.offsets[E][d].get()),
				-1.0
			);
		}
	}

	// ======= S1: collectGlobal* variants for CSR-direct assembly =======

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalZeroNormalIBFIOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			auto op = buildByMult<FES,BF>(
				MInv[f]->SpMat(), buildZeroNormalIBFISubOperator<BF>(f)->SpMat(), fes_);
			for (auto d : { X, Y, Z }) {
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d].get()), -1.0);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalOneNormalIBFIOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue;
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(MInv[f]->SpMat(), buildOneNormalIBFISubOperator<BF>(altField(f), { x })->SpMat(), fes_);
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][y].get(), *globalId.offsets[altField(f)][z].get()),
					1.0 - double(f) * 2.0);
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][z].get(), *globalId.offsets[altField(f)][y].get()),
					-1.0 + double(f) * 2.0);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalTwoNormalIBFIOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				if (d >= dim) continue;
				for (auto d2{ X }; d2 <= Z; d2++) {
					if (d2 >= dim) continue;
					auto op = buildByMult<FES,BF>(MInv[f]->SpMat(), buildTwoNormalIBFISubOperator<BF>(f, { d, d2 })->SpMat(), fes_);
					collectBlockPlacement(op->SpMat(), blocks,
						std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d2].get()), 1.0);
				}
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalDirectionalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue;
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(
					MInv[f]->SpMat(), buildDerivativeSubOperator<BF>(x)->SpMat(), fes_);
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][z].get(), *globalId.offsets[altField(f)][y].get()),
					1.0 - double(f) * 2.0);
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][y].get(), *globalId.offsets[altField(f)][z].get()),
					-1.0 + double(f) * 2.0);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalZeroNormalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			auto op = buildByMult<FES,BF>(
				MInv[f]->SpMat(), buildZeroNormalSubOperator<BF>(f)->SpMat(), fes_);
			for (auto d : { X, Y, Z }) {
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d].get()), -1.0);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalOneNormalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto x{ X }; x <= Z; x++) {
				if (x >= dim) continue;
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult<FES,BF>(
					MInv[f]->SpMat(), buildOneNormalSubOperator<BF>(altField(f), { x })->SpMat(), fes_);
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][y].get(), *globalId.offsets[altField(f)][z].get()),
					1.0 - double(f) * 2.0);
				collectBlockPlacement(op->SpMat(), blocks,
					std::make_pair(*globalId.offsets[f][z].get(), *globalId.offsets[altField(f)][y].get()),
					-1.0 + double(f) * 2.0);
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalTwoNormalOperators(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		const int dim = meshDimension();
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs());
		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				if (d >= dim) continue;
				for (auto d2{ X }; d2 <= Z; d2++) {
					if (d2 >= dim) continue;
					auto op = buildByMult<FES,BF>(
						MInv[f]->SpMat(), buildTwoNormalSubOperator<BF>(f, {d, d2})->SpMat(), fes_);
					collectBlockPlacement(op->SpMat(), blocks,
						std::make_pair(*globalId.offsets[f][d].get(), *globalId.offsets[f][d2].get()), 1.0);
				}
			}
		}
	}

	template <typename FES>
	template <typename BF>
	void DGOperatorFactory<FES>::collectGlobalConductiveOperator(std::vector<CSRBlockPlacement>& blocks, const std::array<std::unique_ptr<BF>, 2>& MInv)
	{
		auto MSig  = buildSigmaMassOperator<BF>();
		auto ASigE = buildByMult<FES, BF>(MInv[E]->SpMat(), MSig->SpMat(), fes_);

		GlobalIndices gid(fes_.GetNDofs(), getAdditionalDofs(), true);
		for (auto d : { X, Y, Z }) {
			collectBlockPlacement(ASigE->SpMat(), blocks,
				std::make_pair(*gid.offsets[E][d].get(), *gid.offsets[E][d].get()), -1.0);
		}
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildSGBCGlobalOperator()
	{
		auto res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * (fes_.GetNDofs() + getAdditionalDofs()));
		auto& interior_marker = pd_.model.getMarker(BdrCond::SGBC, true);
		auto& boundary_marker = pd_.model.getMarker(BdrCond::SGBC, false);

		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
			auto MInv = buildMaxwellInverseMassMatrixOperator<ParBilinearForm>();
			if (interior_marker.Size() != 0 && interior_marker.Sum() != 0) {
				this->template addGlobalSourceFaceIBFIOneNormalOperators<ParBilinearForm>(res.get(), interior_marker, MInv);
				this->template addGlobalSourceFaceIBFIZeroNormalOperators<ParBilinearForm>(res.get(), interior_marker, MInv);
				this->template addGlobalSourceFaceIBFITwoNormalOperators<ParBilinearForm>(res.get(), interior_marker, MInv);
			}
			if (boundary_marker.Size() != 0 && boundary_marker.Sum() != 0) {
				this->template addGlobalBoundarySourceFaceIBFIOneNormalOperators<ParBilinearForm>(res.get(), boundary_marker, MInv);
				this->template addGlobalBoundarySourceFaceIBFIZeroNormalOperators<ParBilinearForm>(res.get(), boundary_marker, MInv);
				this->template addGlobalBoundarySourceFaceIBFITwoNormalOperators<ParBilinearForm>(res.get(), boundary_marker, MInv);
			}
		} else {
			auto MInvSerial = buildMaxwellInverseMassMatrixOperator<BilinearForm>();
			if (interior_marker.Size() != 0 && interior_marker.Sum() != 0) {
				this->template addGlobalSourceFaceIBFIOneNormalOperators<BilinearForm>(res.get(), interior_marker, MInvSerial);
				this->template addGlobalSourceFaceIBFIZeroNormalOperators<BilinearForm>(res.get(), interior_marker, MInvSerial);
				this->template addGlobalSourceFaceIBFITwoNormalOperators<BilinearForm>(res.get(), interior_marker, MInvSerial);
			}
			if (boundary_marker.Size() != 0 && boundary_marker.Sum() != 0) {
				this->template addGlobalBoundarySourceFaceIBFIOneNormalOperators<BilinearForm>(res.get(), boundary_marker, MInvSerial);
				this->template addGlobalBoundarySourceFaceIBFIZeroNormalOperators<BilinearForm>(res.get(), boundary_marker, MInvSerial);
				this->template addGlobalBoundarySourceFaceIBFITwoNormalOperators<BilinearForm>(res.get(), boundary_marker, MInvSerial);
			}
		}

		res->Finalize();
		return res;
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildTFSFGlobalOperator()
	{
		return buildSourceFaceOperator(BdrCond::TotalFieldIn);
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildSourceFaceOperator(BdrCond filter)
	{
		// Look up the marker from the correct model map depending on boundary condition type
		if (filter == BdrCond::TotalFieldIn) {
			auto& tfsfMap = pd_.model.getTotalFieldScatteredFieldToMarker();
			auto it = tfsfMap.find(BdrCond::TotalFieldIn);
			if (it != tfsfMap.end()) {
				return buildSourceFaceOperator(it->second);
			}
		} else if (filter == BdrCond::SGBC) {
			auto& sgbcMap = pd_.model.getSGBCToMarker();
			auto it = sgbcMap.find(BdrCond::SGBC);
			if (it != sgbcMap.end()) {
				return buildSourceFaceOperator(it->second);
			}
		} else {
			auto& intBdrMap = pd_.model.getInteriorBoundaryToMarker();
			auto it = intBdrMap.find(filter);
			if (it != intBdrMap.end()) {
				return buildSourceFaceOperator(it->second);
			}
		}
		// No marker found — return empty finalized matrix
		auto res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * (fes_.GetNDofs() + getAdditionalDofs()));
		res->Finalize();
		return res;
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildSourceFaceOperator(mfem::Array<int>& marker)
	{
		std::unique_ptr<SparseMatrix> res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * (fes_.GetNDofs() + getAdditionalDofs()));

		// S5: Build M^{-1} once and share across all sub-operator assemblies.
		auto MInv = buildMaxwellInverseMassMatrixOperator<ParBilinearForm>();

		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
			// Interior-boundary faces (both elems local): IBFI path.
			this->template addGlobalSourceFaceIBFIOneNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
			this->template addGlobalSourceFaceIBFIZeroNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
			this->template addGlobalSourceFaceIBFITwoNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
			// MPI partition faces (Elem2 on another rank): same Jump source via BFI,
			// mirroring buildSGBCGlobalOperator(). Without this, tagged TFSF faces that
			// appear as mesh boundaries on a rank inject nothing.
			this->template addGlobalBoundarySourceFaceIBFIOneNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
			this->template addGlobalBoundarySourceFaceIBFIZeroNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
			this->template addGlobalBoundarySourceFaceIBFITwoNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
		} else {
			auto MInvSerial = buildMaxwellInverseMassMatrixOperator<BilinearForm>();
			this->template addGlobalSourceFaceIBFIOneNormalOperators<BilinearForm>(res.get(), marker, MInvSerial);
			this->template addGlobalSourceFaceIBFIZeroNormalOperators<BilinearForm>(res.get(), marker, MInvSerial);
			this->template addGlobalSourceFaceIBFITwoNormalOperators<BilinearForm>(res.get(), marker, MInvSerial);
		}

		res->Finalize();
		return res;
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildGlobalOperator()
	{

		if constexpr (std::is_same_v<FES, ParFiniteElementSpace>) {
			fes_.ExchangeFaceNbrData();
    	}

		const int globalRows = 6 * fes_.GetNDofs();
		const int globalCols = 6 * (fes_.GetNDofs() + getAdditionalDofs());

		// S1: Collect all sub-operator block placements, then merge into CSR directly
		// (avoids LIL intermediate representation and its 2x peak memory during Finalize).
		std::vector<CSRBlockPlacement> blocks;

		// S5: Build M^{-1} once and share across all sub-operator assemblies.
		auto MInv = buildMaxwellInverseMassMatrixOperator<ParBilinearForm>();

		if (pd_.model.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

			std::chrono::high_resolution_clock::time_point startTime;
			#ifdef SHOW_TIMER_INFORMATION
			if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
				startTime = std::chrono::high_resolution_clock::now() ;
			}
			#endif

			#ifdef SHOW_TIMER_INFORMATION
			if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
						std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
			}
			#endif

				this->template	collectGlobalOneNormalIBFIOperators<ParBilinearForm>(blocks, MInv);

			#ifdef SHOW_TIMER_INFORMATION
			if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
				startTime = std::chrono::high_resolution_clock::now();
				std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
			}
			#endif

				this->template	collectGlobalZeroNormalIBFIOperators<ParBilinearForm>(blocks, MInv);

			#ifdef SHOW_TIMER_INFORMATION
			if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
				startTime = std::chrono::high_resolution_clock::now();
				std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
			}
			#endif

				this->template	collectGlobalTwoNormalIBFIOperators<ParBilinearForm>(blocks, MInv);

			#ifdef SHOW_TIMER_INFORMATION
			if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			}
			#endif

		}
		else{

			#ifdef SHOW_TIMER_INFORMATION		
					if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
							std::cout << "No Interior Boundary Operators to Assemble." << std::endl;
					}
			#endif

		}

		std::chrono::high_resolution_clock::time_point startTime;
		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
		}
		#endif

		this->template	collectGlobalDirectionalOperators<ParBilinearForm>(blocks, MInv);

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
		}
		#endif

		this->template	collectGlobalOneNormalOperators<ParBilinearForm>(blocks, MInv);

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
		}
		#endif

		this->template	collectGlobalZeroNormalOperators<ParBilinearForm>(blocks, MInv);

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
		}
		#endif

		this->template	collectGlobalTwoNormalOperators<ParBilinearForm>(blocks, MInv);

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Conductivity Operators" << std::endl;
		}
		#endif

		this->template  collectGlobalConductiveOperator<ParBilinearForm>(blocks, MInv);

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "All operators assembled. Merging into global CSR matrix." << std::endl;
		}
		#endif

		auto res = mergeBlocksToCSR(blocks, globalRows, globalCols);

		// Free sub-operator blocks before applying threshold.
		blocks.clear();

		// Threshold set to sqrt(eps_machine) ~ 1e-8, the standard criterion for distinguishing
		// genuine matrix entries from floating-point assembly noise relative to unit-scale quantities.
		auto threshold = 1e-8;
		res->Threshold(threshold);

		if(this->pd_.opts.export_evolution_operator){
			if(Mpi::WorldSize() > 1){
				std::cout << "---------------------------------------------------------------" << std::endl;
				std::cout << "--EXPORTING OPERATOR ONLY CURRENTLY WORKS IN SINGLE RANK SIMS--" << std::endl;
				std::cout << "-----THE SPATIAL EVOLUTION OPERATOR IS NOT BEING EXPORTED------" << std::endl;
				std::cout << "---------------------------------------------------------------" << std::endl;
				return res;
			}
			std::filesystem::path export_dir = std::filesystem::path("Exports") / "Operators" / this->pd_.model.meshName_;

			if (!std::filesystem::exists(export_dir))
			{
				std::filesystem::create_directories(export_dir);
			}

			std::filesystem::path file_path = export_dir / (this->pd_.model.meshName_ + "_global.csr");

			std::ofstream ofs(file_path);
			if (!ofs.is_open())
			{
				throw std::runtime_error("Could not open file for writing: " + file_path.string());
			}

			res->PrintCSR2(ofs);
			ofs.close();

			std::cout << "Global operator exported to " << file_path << std::endl;
		}


		return res;
	}

	namespace {

	/// Bagci SC-PML diagonal tensor entry for field/aux component u.
	enum class SCPMLTensorKind {
		A,       ///< a_uu = κ_v κ_w / κ_u  (Phase 2 LHS; unused in κ≡1 assembly)
		B,       ///< b_uu
		C,       ///< c_uu
		D,       ///< d_uu = σ_u / κ_u
		InvKappa ///< 1/κ_u
	};

	class SCPMLTensorCoefficient : public mfem::Coefficient {
	public:
		SCPMLTensorCoefficient(const PMLProfileData& profiles, Direction comp,
		                       SCPMLTensorKind kind)
			: profiles_(profiles), comp_(comp), kind_(kind)
		{
		}

		double Eval(mfem::ElementTransformation& T,
		            const mfem::IntegrationPoint& ip) override
		{
			double sig[3];
			double kap[3];
			for (int d = 0; d < 3; ++d) {
				PMLDirectionProfiles out;
				profiles_.evaluateAtTransform(T, ip, static_cast<Direction>(d), out);
				sig[d] = out.sigma;
				kap[d] = std::max(out.kappa, 1e-30);
			}
			const int u = static_cast<int>(comp_);
			const int v = (u + 1) % 3;
			const int w = (u + 2) % 3;
			const double a = kap[v] * kap[w] / kap[u];
			const double b =
				(sig[v] * kap[w] + sig[w] * kap[v] - a * sig[u]) / kap[u];
			const double c = sig[v] * sig[w] - b * sig[u];
			const double d = sig[u] / kap[u];
			switch (kind_) {
			case SCPMLTensorKind::A:
				return a;
			case SCPMLTensorKind::B:
				return b;
			case SCPMLTensorKind::C:
				return c;
			case SCPMLTensorKind::D:
				return d;
			case SCPMLTensorKind::InvKappa:
				return 1.0 / kap[u];
			}
			return 0.0;
		}

	private:
		const PMLProfileData& profiles_;
		Direction comp_;
		SCPMLTensorKind kind_;
	};

	} // namespace

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildMarkedMassOperator(
		mfem::Coefficient& coeff, mfem::Array<int>& attr_marker)
	{
		auto bf = std::make_unique<BF>(&fes_);
		bf->AddDomainIntegrator(new MassIntegrator(coeff), attr_marker);
		bf->Assemble();
		bf->Finalize();
		return bf;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildMarkedInverseMassOperator(
		mfem::Coefficient& coeff, mfem::Array<int>& attr_marker)
	{
		auto bf = std::make_unique<BF>(&fes_);
		bf->AddDomainIntegrator(
			new InverseIntegrator(new MassIntegrator(coeff)), attr_marker);
		bf->Assemble();
		bf->Finalize();
		return bf;
	}

	template <typename FES>
	void DGOperatorFactory<FES>::buildSCPMLOperators(
		const SCPMLLayout& layout,
		std::unique_ptr<mfem::SparseMatrix>& ade_operator,
		std::array<std::unique_ptr<mfem::SparseMatrix>, 3>& curl_delta)
	{
		ade_operator.reset();
		for (auto& d : curl_delta) {
			d.reset();
		}
		if (layout.nAux() == 0) {
			return;
		}
		const PMLProfileData* profiles = pd_.model.getPMLProfileData();
		if (!profiles) {
			throw std::runtime_error(
				"buildSCPMLOperators requires initialized PML profile data.");
		}

		bool needs_a = false;
		for (const auto& props : pd_.model.getPMLProperties()) {
			if (props.kappa_max > 1.0 + 1e-12) {
				needs_a = true;
				break;
			}
		}

		const int ndofs = fes_.GetNDofs();
		const int n_aux = layout.nAux();
		const int globalRows = 6 * ndofs + n_aux;
		const int globalCols = globalRows;

		mfem::Array<int> pml_marker = pd_.model.buildPMLVolumeMarker();
		auto MInv = buildMaxwellInverseMassMatrixOperator<ParBilinearForm>();

		mfem::ConstantCoefficient one(1.0);
		auto Munit = buildMarkedMassOperator<ParBilinearForm>(one, pml_marker);
		auto S_unit_E = buildByMult<FES, ParBilinearForm>(
			MInv[E]->SpMat(), Munit->SpMat(), fes_);
		auto S_unit_H = buildByMult<FES, ParBilinearForm>(
			MInv[H]->SpMat(), Munit->SpMat(), fes_);

		std::vector<CSRBlockPlacement> blocks;

		for (Direction u = X; u <= Z; ++u) {
			SCPMLTensorCoefficient c_a(*profiles, u, SCPMLTensorKind::A);
			SCPMLTensorCoefficient c_b(*profiles, u, SCPMLTensorKind::B);
			SCPMLTensorCoefficient c_c(*profiles, u, SCPMLTensorKind::C);
			SCPMLTensorCoefficient c_d(*profiles, u, SCPMLTensorKind::D);
			SCPMLTensorCoefficient c_invk(*profiles, u, SCPMLTensorKind::InvKappa);

			auto Mb = buildMarkedMassOperator<ParBilinearForm>(c_b, pml_marker);
			auto Mc = buildMarkedMassOperator<ParBilinearForm>(c_c, pml_marker);
			auto Md = buildMarkedMassOperator<ParBilinearForm>(c_d, pml_marker);
			auto Minvk = buildMarkedMassOperator<ParBilinearForm>(c_invk, pml_marker);

			std::unique_ptr<ParBilinearForm> A_b_E;
			std::unique_ptr<ParBilinearForm> A_b_H;
			std::unique_ptr<ParBilinearForm> A_c_E;
			std::unique_ptr<ParBilinearForm> A_c_H;
			if (needs_a) {
				auto MaInv = buildMarkedInverseMassOperator<ParBilinearForm>(
					c_a, pml_marker);
				A_b_E = buildByMult<FES, ParBilinearForm>(
					MaInv->SpMat(), Mb->SpMat(), fes_);
				A_b_H = buildByMult<FES, ParBilinearForm>(
					MaInv->SpMat(), Mb->SpMat(), fes_);
				A_c_E = buildByMult<FES, ParBilinearForm>(
					MaInv->SpMat(), Mc->SpMat(), fes_);
				A_c_H = buildByMult<FES, ParBilinearForm>(
					MaInv->SpMat(), Mc->SpMat(), fes_);

				// Delta = MaInv*M − MInv*M so out += Delta*out ⇒ MaInv*M*out on PML.
				auto R = buildByMult<FES, ParBilinearForm>(
					MaInv->SpMat(), Munit->SpMat(), fes_);
				auto delta = std::make_unique<SparseMatrix>(R->SpMat());
				delta->Add(-1.0, S_unit_E->SpMat());
				delta->Finalize();
				if (delta->NumNonZeroElems() > 0) {
					curl_delta[u] = std::move(delta);
				}
			} else {
				A_b_E = buildByMult<FES, ParBilinearForm>(
					MInv[E]->SpMat(), Mb->SpMat(), fes_);
				A_b_H = buildByMult<FES, ParBilinearForm>(
					MInv[H]->SpMat(), Mb->SpMat(), fes_);
				A_c_E = buildByMult<FES, ParBilinearForm>(
					MInv[E]->SpMat(), Mc->SpMat(), fes_);
				A_c_H = buildByMult<FES, ParBilinearForm>(
					MInv[H]->SpMat(), Mc->SpMat(), fes_);
			}

			// P rows always use unit Maxwell MInv (no a on ∂t P).
			auto A_d_E = buildByMult<FES, ParBilinearForm>(
				MInv[E]->SpMat(), Md->SpMat(), fes_);
			auto A_d_H = buildByMult<FES, ParBilinearForm>(
				MInv[H]->SpMat(), Md->SpMat(), fes_);
			auto A_invk_E = buildByMult<FES, ParBilinearForm>(
				MInv[E]->SpMat(), Minvk->SpMat(), fes_);
			auto A_invk_H = buildByMult<FES, ParBilinearForm>(
				MInv[H]->SpMat(), Minvk->SpMat(), fes_);

			const int pe = layout.pEOffset(u);
			const int ph = layout.pHOffset(u);
			const int e_off = u * ndofs;
			const int h_off = (3 + u) * ndofs;

			collectBlockPlacement(A_b_E->SpMat(), blocks, e_off, e_off, -1.0);
			collectBlockPlacement(A_c_E->SpMat(), blocks, e_off, pe, -1.0);
			collectBlockPlacement(A_b_H->SpMat(), blocks, h_off, h_off, -1.0);
			collectBlockPlacement(A_c_H->SpMat(), blocks, h_off, ph, -1.0);

			collectBlockPlacement(A_invk_E->SpMat(), blocks, pe, e_off, 1.0);
			collectBlockPlacement(A_d_E->SpMat(), blocks, pe, pe, -1.0);
			collectBlockPlacement(A_invk_H->SpMat(), blocks, ph, h_off, 1.0);
			collectBlockPlacement(A_d_H->SpMat(), blocks, ph, ph, -1.0);
		}

		(void)S_unit_H;

		int sum_nnz = 0;
		for (const auto& bp : blocks) {
			if (bp.block) {
				sum_nnz += bp.block->NumNonZeroElems();
			}
		}
		if (sum_nnz <= 0) {
			blocks.clear();
			for (auto& d : curl_delta) {
				d.reset();
			}
			std::cout << "[PML] Rank " << Mpi::WorldRank()
			          << ": no local PML volume — SC-PML ADE operator omitted"
			          << std::endl;
			return;
		}

		ade_operator = mergeBlocksToCSR(blocks, globalRows, globalCols);
		blocks.clear();
		ade_operator->Threshold(1e-8);

		std::cout << "[PML] Rank " << Mpi::WorldRank()
		          << ": SC-PML ADE operator " << globalRows << " x " << globalCols
		          << ", nnz=" << ade_operator->NumNonZeroElems()
		          << (needs_a ? " (MaInv damping + curl a-rescale)"
		                      : " (κ≡1 unit MInv)")
		          << std::endl;
	}

} // namespace maxwell
