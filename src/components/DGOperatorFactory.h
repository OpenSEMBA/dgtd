#pragma once

#include "ProblemDescription.h"
#include "PMLAuxLayout.h"
#include "PMLCoefficients.h"
#include "PMLDGHelpers.h"
#include "PMLSignTest.h"
#include "PMLOperatorAudit.h"

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
			throw std::runtime_error(
				"mergeBlocksToCSR: total NNZ (" + std::to_string(totalNNZ_64) +
				") is out of range for SparseMatrix indexing.");
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
		std::unique_ptr<mfem::SparseMatrix> buildPMLOperator(const PMLAuxLayout& layout);

	private:
		template <typename BF>
		std::unique_ptr<BF> buildPMLDomainMassSubOperator(
			mfem::Coefficient& coeff, mfem::Array<int>& pml_marker, int ir_order);

		template <typename BF>
		std::unique_ptr<BF> buildPMLDomainDerivativeSubOperator(
			mfem::Coefficient& coeff, Direction deriv_dir,
			mfem::Array<int>& pml_marker);

		template <typename BF>
		std::unique_ptr<BF> buildPMLDomainOneNormalSubOperator(
			mfem::Coefficient& coeff, const FieldType& field,
			const std::vector<Direction>& dir_terms, mfem::Array<int>& pml_marker);

		template <typename BF>
		std::unique_ptr<BF> buildPMLDomainZeroNormalSubOperator(
			mfem::Coefficient& coeff, mfem::Array<int>& pml_marker);

		template <typename BF>
		std::unique_ptr<BF> buildPMLDomainTwoNormalSubOperator(
			mfem::Coefficient& coeff, const FieldType& field,
			const std::vector<Direction>& dir_terms, mfem::Array<int>& pml_marker);

		template <typename BF>
		std::unique_ptr<BF> buildPMLScalarInverseMassSubOperator(
			mfem::Array<int>& pml_marker, int ir_order);

		struct PMLOperatorAssemblyStats {
			int driver_volume_nnz = 0;
			int driver_face_nnz = 0;
			int driver_upwind_nnz = 0;
			int correction_nnz = 0;
		};

		void collectPMLUpwindDriverBlocks(
			std::vector<CSRBlockPlacement>& blocks,
			int psi_row,
			FieldType in_field,
			Direction in_c,
			const mfem::SparseMatrix& zero_op,
			const std::vector<std::pair<Direction, Direction>>& two_dir_pairs,
			const std::vector<const mfem::SparseMatrix*>& two_ops,
			const GlobalIndices& globalId,
			PMLOperatorAssemblyStats& stats);

		void collectPMLComponentDriverBlocks(
			std::vector<CSRBlockPlacement>& blocks,
			int psi_row,
			FieldType in_field,
			Direction in_c,
			Direction stretch_dir,
			double weight,
			const mfem::SparseMatrix& vol_op,
			const mfem::SparseMatrix& face_op,
			const GlobalIndices& globalId,
			PMLOperatorAssemblyStats& stats);

		void collectPMLOperatorBlocks(
			std::vector<CSRBlockPlacement>& blocks,
			const PMLAuxLayout& layout,
			const PMLProfileData& profiles,
			mfem::Array<int>& pml_marker,
			const std::unique_ptr<ParBilinearForm>& MInvScalar,
			const std::array<std::unique_ptr<ParBilinearForm>, 2>& MInvMaxwell,
			PMLOperatorAssemblyStats& stats);

		void auditGlobalOperatorCurl(const mfem::SparseMatrix& global_op);

		void auditPMLOperatorCurl(const PMLAuxLayout& layout,
		                          const mfem::SparseMatrix& pml_op,
		                          mfem::Array<int>& pml_marker,
		                          const PMLProfileData& profiles);

		ProblemDescription pd_;
		FES fes_;

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
			if (kv.first == BdrCond::SGBC || kv.first == BdrCond::PML_NONE)
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
			if (kv.first == BdrCond::SGBC || kv.first == BdrCond::PML_NONE)
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
			if (kv.first == BdrCond::SGBC || kv.first == BdrCond::PML_NONE)
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
		if (isPMLOperatorAuditEnabled() && Mpi::WorldRank() == 0) {
			pmlAuditLog("GlobalDirectional placements=" + std::to_string(blocks.size()));
		}

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
		}
		#endif

		this->template	collectGlobalOneNormalOperators<ParBilinearForm>(blocks, MInv);
		if (isPMLOperatorAuditEnabled() && Mpi::WorldRank() == 0) {
			pmlAuditLog("GlobalOneNormal placements=" + std::to_string(blocks.size()));
		}

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
		}
		#endif

		this->template	collectGlobalZeroNormalOperators<ParBilinearForm>(blocks, MInv);
		if (isPMLOperatorAuditEnabled() && Mpi::WorldRank() == 0) {
			pmlAuditLog("GlobalZeroNormal placements=" + std::to_string(blocks.size()) +
			            " upwind_alpha=" + std::to_string(pd_.opts.alpha));
		}

		#ifdef SHOW_TIMER_INFORMATION
		if (!pd_.opts.is_sgbc_solver && Mpi::WorldRank() == 0){
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			startTime = std::chrono::high_resolution_clock::now();
			std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
		}
		#endif

		this->template	collectGlobalTwoNormalOperators<ParBilinearForm>(blocks, MInv);
		if (isPMLOperatorAuditEnabled() && Mpi::WorldRank() == 0) {
			pmlAuditLog("GlobalTwoNormal placements=" + std::to_string(blocks.size()));
		}

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

		if (isPMLOperatorAuditEnabled()) {
			auditGlobalOperatorCurl(*res);
		}

		return res;
	}

	template <typename FES>
	void DGOperatorFactory<FES>::auditGlobalOperatorCurl(
		const mfem::SparseMatrix& global_op)
	{
		if (Mpi::WorldRank() != 0) {
			return;
		}
		printIntegratorBaseline1DTE();

		const int ndofs = fes_.GetNDofs();
		const int nbrDofs = getAdditionalDofs();
		const int blockSize = ndofs + nbrDofs;
		const int dim = meshDimension();
		GlobalIndices gid(ndofs, nbrDofs, true);

		const int hz_row = gid.offsets[H][Z]->rowStartOffset;
		const int ey_col = gid.offsets[E][Y]->colStartOffset;

		const double fn_hz_ey_full = frobeniusBlockNorm(
			global_op, hz_row, hz_row + ndofs, ey_col, ey_col + blockSize);

		pmlAuditLog("global full nnz=" + std::to_string(global_op.NumNonZeroElems()) +
		            " ||Hz<-Ey||_F=" + std::to_string(fn_hz_ey_full));

		if (dim < 1 || !pd_.model.hasPML()) {
			return;
		}

		mfem::Array<int> pml_marker = pd_.model.buildPMLVolumeMarker();
		std::vector<int> pml_dofs;
		mfem::Mesh& mesh = *fes_.GetMesh();
		const mfem::Table& elem_dof_table = fes_.GetElementToDofTable();
		for (int el = 0; el < mesh.GetNE(); ++el) {
			const int attr = mesh.GetAttribute(el);
			if (attr < 1 || attr > pml_marker.Size() || pml_marker[attr - 1] == 0) {
				continue;
			}
			const int* dofs = elem_dof_table.GetRow(el);
			for (int i = 0; i < elem_dof_table.RowSize(el); ++i) {
				pml_dofs.push_back(dofs[i]);
			}
		}
		std::sort(pml_dofs.begin(), pml_dofs.end());
		pml_dofs.erase(std::unique(pml_dofs.begin(), pml_dofs.end()), pml_dofs.end());
		if (!pml_dofs.empty()) {
			const double fn_hz_ey_pml = frobeniusBlockNormOnRows(
				global_op, pml_dofs, hz_row, ey_col, ey_col + blockSize);
			pmlAuditLog("||Hz<-Ey||_F on PML element DOFs=" + std::to_string(fn_hz_ey_pml));
		}

		auto MInv = buildMaxwellInverseMassMatrixOperator<ParBilinearForm>();
		const Direction x = X;
		const Direction y = static_cast<Direction>((x + 1) % 3);
		const Direction z = static_cast<Direction>((x + 2) % 3);

		auto op_dir = buildByMult<FES, ParBilinearForm>(
			MInv[H]->SpMat(), buildDerivativeSubOperator<ParBilinearForm>(x)->SpMat(), fes_);
		auto op_one = buildByMult<FES, ParBilinearForm>(
			MInv[H]->SpMat(), buildOneNormalSubOperator<ParBilinearForm>(E, { x })->SpMat(), fes_);

		std::vector<CSRBlockPlacement> dx_blocks;
		collectBlockPlacement(op_dir->SpMat(), dx_blocks,
			std::make_pair(*gid.offsets[H][z].get(), *gid.offsets[E][y].get()), -1.0);
		collectBlockPlacement(op_one->SpMat(), dx_blocks,
			std::make_pair(*gid.offsets[H][z].get(), *gid.offsets[E][y].get()), -1.0);

		auto dx_merged = mergeBlocksToCSR(dx_blocks, 6 * ndofs, 6 * blockSize);
		const double fn_dx = frobeniusBlockNorm(
			*dx_merged, hz_row, hz_row + ndofs, ey_col, ey_col + blockSize);

		pmlAuditLog("global_Dx+OneNormal_x nnz=" + std::to_string(dx_merged->NumNonZeroElems()) +
		            " ||Hz<-Ey||_F=" + std::to_string(fn_dx) +
		            " ratio_vs_full=" + std::to_string(fn_dx / std::max(fn_hz_ey_full, 1e-30)));

		if (shouldExportPMLOperatorAuditCSR() && Mpi::WorldSize() == 1) {
			const std::string audit_subdir = pd_.model.meshName_ + "_audit";
			const std::filesystem::path export_dir =
				std::filesystem::path("Exports") / "Operators" / audit_subdir;
			std::filesystem::create_directories(export_dir);
			exportAuditCSR(global_op, export_dir.string(), "global_full.csr");
			exportAuditCSR(*dx_merged, export_dir.string(), "global_Dx_only.csr");
		}
	}

	template <typename FES>
	void DGOperatorFactory<FES>::auditPMLOperatorCurl(
		const PMLAuxLayout& layout,
		const mfem::SparseMatrix& pml_op,
		mfem::Array<int>& pml_marker,
		const PMLProfileData& profiles)
	{
		if (Mpi::WorldRank() != 0 || layout.numStretchDirections() == 0) {
			return;
		}

		const int ndofs = fes_.GetNDofs();
		const int nbrDofs = getAdditionalDofs();
		const int blockSize = ndofs + nbrDofs;
		const Direction stretch_dir = layout.stretchDirection(0);
		const int dim = meshDimension();
		if (stretch_dir < 0 || stretch_dir >= dim) {
			return;
		}

		PMLProfileCoefficient sigma_coeff(profiles, stretch_dir, PMLProfileCoefficient::Kind::Sigma);
		auto MInvScalar = buildPMLScalarInverseMassSubOperator<ParBilinearForm>(
			pml_marker, profiles.feOrder() + 1);

		const std::vector<Direction> dir_terms = { stretch_dir };
		auto psi_e_volume = buildPMLDomainDerivativeSubOperator<ParBilinearForm>(
			sigma_coeff, stretch_dir, pml_marker);
		auto psi_e_face = buildPMLDomainOneNormalSubOperator<ParBilinearForm>(
			sigma_coeff, altField(H), dir_terms, pml_marker);
		auto psi_e_vol_op = buildByMult<FES, ParBilinearForm>(
			MInvScalar->SpMat(), psi_e_volume->SpMat(), fes_);
		auto psi_e_face_op = buildByMult<FES, ParBilinearForm>(
			MInvScalar->SpMat(), psi_e_face->SpMat(), fes_);

		Direction in_c = X;
		double weight = 0.0;
		const Direction h_comp = Z;
		if (!pmlPsiEDriverCoupling(h_comp, stretch_dir, in_c, weight)) {
			pmlAuditLog("psi^E driver coupling inactive for 1D TE audit");
			return;
		}

		GlobalIndices globalId(ndofs, nbrDofs, true);
		const int psi_row = layout.psiEOffset(stretch_dir, h_comp);
		const int ey_col = globalId.offsets[E][in_c]->colStartOffset;
		const int hz_row = globalId.offsets[H][h_comp]->rowStartOffset;
		const int psi_col = 6 * blockSize + (psi_row - 6 * ndofs);

		std::vector<CSRBlockPlacement> driver_blocks;
		collectBlockPlacement(
			psi_e_vol_op->SpMat(), driver_blocks, psi_row, ey_col, weight);
		collectBlockPlacement(
			psi_e_face_op->SpMat(), driver_blocks, psi_row, ey_col, -weight);

		auto driver_merged = mergeBlocksToCSR(
			driver_blocks, 6 * ndofs + layout.nAux(), 6 * blockSize + layout.nAux());

		const double fn_psi_ey = frobeniusBlockNorm(
			*driver_merged, psi_row, psi_row + ndofs, ey_col, ey_col + blockSize);
		const double fn_hz_psi = frobeniusBlockNorm(
			pml_op, hz_row, hz_row + ndofs, psi_col, psi_col + ndofs);
		const double fn_psi_psi = frobeniusBlockNorm(
			pml_op, psi_row, psi_row + ndofs, psi_col, psi_col + ndofs);

		pmlAuditLog("PML driver vol nnz=" + std::to_string(psi_e_vol_op->SpMat().NumNonZeroElems()) +
		            " face nnz=" + std::to_string(psi_e_face_op->SpMat().NumNonZeroElems()) +
		            " w=" + std::to_string(weight));
		pmlAuditLog("||psi<-Ey||_F=" + std::to_string(fn_psi_ey) +
		            " ||Hz<-psi||_F=" + std::to_string(fn_hz_psi) +
		            " ||psi<-psi_mass||_F=" + std::to_string(fn_psi_psi));

		std::vector<int> pml_dofs;
		mfem::Mesh& mesh = *fes_.GetMesh();
		const mfem::Table& elem_dof_table = fes_.GetElementToDofTable();
		for (int el = 0; el < mesh.GetNE(); ++el) {
			const int attr = mesh.GetAttribute(el);
			if (attr < 1 || attr > pml_marker.Size() || pml_marker[attr - 1] == 0) {
				continue;
			}
			const int* dofs = elem_dof_table.GetRow(el);
			for (int i = 0; i < elem_dof_table.RowSize(el); ++i) {
				pml_dofs.push_back(dofs[i]);
			}
		}
		std::sort(pml_dofs.begin(), pml_dofs.end());
		pml_dofs.erase(std::unique(pml_dofs.begin(), pml_dofs.end()), pml_dofs.end());

		if (!pml_dofs.empty()) {
			pmlAuditLog("overlap note: compare global ||Hz<-Ey|| on PML DOFs vs ||psi<-Ey|| "
			            "(duplicate curl if same order)");
		}

		if (shouldExportPMLOperatorAuditCSR() && Mpi::WorldSize() == 1) {
			const std::string audit_subdir = pd_.model.meshName_ + "_audit";
			const std::filesystem::path export_dir =
				std::filesystem::path("Exports") / "Operators" / audit_subdir;
			std::filesystem::create_directories(export_dir);
			exportAuditCSR(*driver_merged, export_dir.string(), "pml_driver_Ey_to_psi.csr");
			exportAuditCSR(pml_op, export_dir.string(), "pml_full.csr");
		}
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildPMLScalarInverseMassSubOperator(
		mfem::Array<int>& pml_marker, int ir_order)
	{
		(void)pml_marker;
		(void)ir_order;
		mfem::ConstantCoefficient one(1.0);
		auto res = std::make_unique<BF>(&fes_);
		res->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildPMLDomainMassSubOperator(
		mfem::Coefficient& coeff, mfem::Array<int>& pml_marker, int ir_order)
	{
		auto res = std::make_unique<BF>(&fes_);
		auto* integ = new MassIntegrator(coeff);
		if (fes_.GetMesh()->GetNE() > 0) {
			const auto geom = fes_.GetMesh()->GetElementGeometry(0);
			integ->SetIntRule(&IntRules.Get(geom, ir_order));
		}
		res->AddDomainIntegrator(integ, pml_marker);
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildPMLDomainDerivativeSubOperator(
		mfem::Coefficient& coeff, Direction deriv_dir,
		mfem::Array<int>& pml_marker)
	{
		auto res = std::make_unique<BF>(&fes_);
		if (deriv_dir >= fes_.GetMesh()->Dimension()) {
			res->Assemble();
			res->Finalize();
			return res;
		}

		auto* integ = new DerivativeIntegrator(coeff, deriv_dir);
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

		res->AddDomainIntegrator(integ, pml_marker);
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildPMLDomainOneNormalSubOperator(
		mfem::Coefficient& coeff, const FieldType& field,
		const std::vector<Direction>& dir_terms, mfem::Array<int>& pml_marker)
	{
		(void)field;
		auto res = std::make_unique<BF>(&fes_);
		res->AddInteriorFaceIntegrator(
			new MaxwellDGCoefficientOneNormalJumpIntegrator(dir_terms, coeff), pml_marker);
		// Optional terminating-boundary faces (S8 via env PML_SIGN_TEST=8 — default off; trial worsened t=20).
		if (getPMLSignTestMode() == PMLSignTestMode::IncludeOuterBdyFace) {
			for (auto& kv : pd_.model.getBoundaryToMarker()) {
				if (kv.first == BdrCond::SGBC) {
					continue;
				}
				res->AddBdrFaceIntegrator(
					new MaxwellDGCoefficientOneNormalJumpIntegrator(dir_terms, coeff),
					kv.second);
			}
		}
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildPMLDomainZeroNormalSubOperator(
		mfem::Coefficient& coeff, mfem::Array<int>& pml_marker)
	{
		auto res = std::make_unique<BF>(&fes_);
		res->AddInteriorFaceIntegrator(
			new MaxwellDGCoefficientZeroNormalJumpIntegrator(coeff), pml_marker);
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildPMLDomainTwoNormalSubOperator(
		mfem::Coefficient& coeff, const FieldType& field,
		const std::vector<Direction>& dir_terms, mfem::Array<int>& pml_marker)
	{
		(void)field;
		auto res = std::make_unique<BF>(&fes_);
		res->AddInteriorFaceIntegrator(
			new MaxwellDGCoefficientTwoNormalJumpIntegrator(dir_terms, coeff), pml_marker);
		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	void DGOperatorFactory<FES>::collectPMLUpwindDriverBlocks(
		std::vector<CSRBlockPlacement>& blocks,
		int psi_row,
		FieldType in_field,
		Direction in_c,
		const mfem::SparseMatrix& zero_op,
		const std::vector<std::pair<Direction, Direction>>& two_dir_pairs,
		const std::vector<const mfem::SparseMatrix*>& two_ops,
		const GlobalIndices& globalId,
		PMLOperatorAssemblyStats& stats)
	{
		collectBlockPlacement(
			zero_op, blocks, psi_row,
			globalId.offsets[in_field][in_c]->colStartOffset, -1.0);
		stats.driver_upwind_nnz += zero_op.NumNonZeroElems();

		for (size_t i = 0; i < two_dir_pairs.size(); ++i) {
			if (two_dir_pairs[i].first != in_c) {
				continue;
			}
			const Direction d2 = two_dir_pairs[i].second;
			collectBlockPlacement(
				*two_ops[i], blocks, psi_row,
				globalId.offsets[in_field][d2]->colStartOffset, 1.0);
			stats.driver_upwind_nnz += two_ops[i]->NumNonZeroElems();
		}
	}

	template <typename FES>
	void DGOperatorFactory<FES>::collectPMLComponentDriverBlocks(
		std::vector<CSRBlockPlacement>& blocks,
		int psi_row,
		FieldType in_field,
		Direction in_c,
		Direction stretch_dir,
		double weight,
		const mfem::SparseMatrix& vol_op,
		const mfem::SparseMatrix& face_op,
		const GlobalIndices& globalId,
		PMLOperatorAssemblyStats& stats)
	{
		(void)stats;
		const PMLSignTestMode mode = getPMLSignTestMode();
		double w = weight;
		if (mode == PMLSignTestMode::NegateDriverWeight) {
			w = -w;
		}

		collectBlockPlacement(
			vol_op, blocks, psi_row,
			globalId.offsets[in_field][in_c]->colStartOffset, w);

		if (mode == PMLSignTestMode::FaceSameAsVolume) {
			collectBlockPlacement(
				face_op, blocks, psi_row,
				globalId.offsets[in_field][in_c]->colStartOffset, w);
		} else if (mode == PMLSignTestMode::FaceCrossColumn) {
			const Direction y = static_cast<Direction>((stretch_dir + 1) % 3);
			const Direction z = static_cast<Direction>((stretch_dir + 2) % 3);
			const FieldType f = in_field;
			collectBlockPlacement(
				face_op, blocks, psi_row,
				globalId.offsets[altField(f)][z]->colStartOffset,
				w * (1.0 - double(f) * 2.0));
			collectBlockPlacement(
				face_op, blocks, psi_row,
				globalId.offsets[altField(f)][y]->colStartOffset,
				w * (-1.0 + double(f) * 2.0));
		} else {
			// Default SBP: face jump opposes volume on the same field column.
			collectBlockPlacement(
				face_op, blocks, psi_row,
				globalId.offsets[in_field][in_c]->colStartOffset, -w);
		}
	}

	template <typename FES>
	void DGOperatorFactory<FES>::collectPMLOperatorBlocks(
		std::vector<CSRBlockPlacement>& blocks,
		const PMLAuxLayout& layout,
		const PMLProfileData& profiles,
		mfem::Array<int>& pml_marker,
		const std::unique_ptr<ParBilinearForm>& MInvScalar,
		const std::array<std::unique_ptr<ParBilinearForm>, 2>& MInvMaxwell,
		PMLOperatorAssemblyStats& stats)
	{
		const int ndofs = fes_.GetNDofs();
		const int nbrDofs = getAdditionalDofs();
		const int blockSize = ndofs + nbrDofs;
		const int dim = meshDimension();
		const int ir_order = profiles.feOrder() + 1;
		GlobalIndices globalId(ndofs, nbrDofs, true);
		const PMLSignTestMode sign_mode = getPMLSignTestMode();
		const double psi_mass_sign =
			(sign_mode == PMLSignTestMode::FlipPsiMassSign) ? 1.0 : -1.0;
		const double h_corr_sign =
			(sign_mode == PMLSignTestMode::FlipCorrections) ? 1.0 : -1.0;
		const double e_corr_sign =
			(sign_mode == PMLSignTestMode::FlipCorrections) ? -1.0 : 1.0;
		const bool pml_upwind = (pd_.opts.alpha > 0.0);

		std::vector<std::pair<Direction, Direction>> two_dir_pairs;
		if (pml_upwind) {
			for (auto d : { X, Y, Z }) {
				if (d >= dim) {
					continue;
				}
				for (auto d2 : { X, Y, Z }) {
					if (d2 >= dim) {
						continue;
					}
					two_dir_pairs.emplace_back(d, d2);
				}
			}
		}

		for (int slot = 0; slot < layout.numStretchDirections(); ++slot) {
			const Direction stretch_dir = layout.stretchDirection(slot);
			if (stretch_dir < 0 || stretch_dir >= dim) {
				continue;
			}

			PMLProfileCoefficient alpha_coeff(profiles, stretch_dir, PMLProfileCoefficient::Kind::Alpha);
			PMLProfileCoefficient sigma_coeff(profiles, stretch_dir, PMLProfileCoefficient::Kind::Sigma);
			PMLProfileCoefficient inv_kappa_coeff(profiles, stretch_dir, PMLProfileCoefficient::Kind::InvKappa);

			const std::vector<Direction> dir_terms = {stretch_dir};

			auto psi_e_volume = buildPMLDomainDerivativeSubOperator<ParBilinearForm>(
				sigma_coeff, stretch_dir, pml_marker);
			auto psi_e_face = buildPMLDomainOneNormalSubOperator<ParBilinearForm>(
				sigma_coeff, altField(H), dir_terms, pml_marker);
			auto psi_h_volume = buildPMLDomainDerivativeSubOperator<ParBilinearForm>(
				sigma_coeff, stretch_dir, pml_marker);
			auto psi_h_face = buildPMLDomainOneNormalSubOperator<ParBilinearForm>(
				sigma_coeff, altField(E), dir_terms, pml_marker);

			auto psi_e_vol_op = buildByMult<FES, ParBilinearForm>(
				MInvScalar->SpMat(), psi_e_volume->SpMat(), fes_);
			auto psi_e_face_op = buildByMult<FES, ParBilinearForm>(
				MInvScalar->SpMat(), psi_e_face->SpMat(), fes_);
			auto psi_h_vol_op = buildByMult<FES, ParBilinearForm>(
				MInvScalar->SpMat(), psi_h_volume->SpMat(), fes_);
			auto psi_h_face_op = buildByMult<FES, ParBilinearForm>(
				MInvScalar->SpMat(), psi_h_face->SpMat(), fes_);

			auto psi_mass = buildPMLDomainMassSubOperator<ParBilinearForm>(
				alpha_coeff, pml_marker, ir_order);
			auto psi_mass_op = buildByMult<FES, ParBilinearForm>(
				MInvScalar->SpMat(), psi_mass->SpMat(), fes_);

			auto field_corr = buildPMLDomainMassSubOperator<ParBilinearForm>(
				inv_kappa_coeff, pml_marker, ir_order);
			auto e_corr_op = buildByMult<FES, ParBilinearForm>(
				MInvMaxwell[E]->SpMat(), field_corr->SpMat(), fes_);
			auto h_corr_op = buildByMult<FES, ParBilinearForm>(
				MInvMaxwell[H]->SpMat(), field_corr->SpMat(), fes_);

			std::unique_ptr<ParBilinearForm> psi_zero_op;
			std::vector<std::unique_ptr<ParBilinearForm>> psi_two_ops_E;
			std::vector<std::unique_ptr<ParBilinearForm>> psi_two_ops_H;
			std::vector<const mfem::SparseMatrix*> psi_two_ptr_E;
			std::vector<const mfem::SparseMatrix*> psi_two_ptr_H;
			if (pml_upwind) {
				ConstantCoefficient upwind_scale(pd_.opts.alpha);
				ProductCoefficient sigma_upwind(sigma_coeff, upwind_scale);
				auto psi_zero_bf = buildPMLDomainZeroNormalSubOperator<ParBilinearForm>(
					sigma_upwind, pml_marker);
				psi_zero_op = buildByMult<FES, ParBilinearForm>(
					MInvScalar->SpMat(), psi_zero_bf->SpMat(), fes_);
				for (const auto& pair : two_dir_pairs) {
					auto two_e_bf = buildPMLDomainTwoNormalSubOperator<ParBilinearForm>(
						sigma_upwind, E, { pair.first, pair.second }, pml_marker);
					auto two_h_bf = buildPMLDomainTwoNormalSubOperator<ParBilinearForm>(
						sigma_upwind, H, { pair.first, pair.second }, pml_marker);
					psi_two_ops_E.push_back(buildByMult<FES, ParBilinearForm>(
						MInvScalar->SpMat(), two_e_bf->SpMat(), fes_));
					psi_two_ops_H.push_back(buildByMult<FES, ParBilinearForm>(
						MInvScalar->SpMat(), two_h_bf->SpMat(), fes_));
					psi_two_ptr_E.push_back(&psi_two_ops_E.back()->SpMat());
					psi_two_ptr_H.push_back(&psi_two_ops_H.back()->SpMat());
				}
			}

			for (Direction h_comp = X; h_comp <= Z; ++h_comp) {
				if (!pmlPsiEComponentActive(h_comp, stretch_dir)) {
					continue;
				}
				Direction in_c = X;
				double weight = 0.0;
				if (!pmlPsiEDriverCoupling(h_comp, stretch_dir, in_c, weight)) {
					continue;
				}

				const int psi_e_row = layout.psiEOffset(stretch_dir, h_comp);
				const int psi_e_col = 6 * blockSize + (psi_e_row - 6 * ndofs);

				collectBlockPlacement(
					psi_mass_op->SpMat(), blocks, psi_e_row, psi_e_col, psi_mass_sign);
				collectPMLComponentDriverBlocks(
					blocks, psi_e_row, E, in_c, stretch_dir, weight,
					psi_e_vol_op->SpMat(), psi_e_face_op->SpMat(),
					globalId, stats);
				stats.driver_volume_nnz += psi_e_vol_op->SpMat().NumNonZeroElems();
				stats.driver_face_nnz += psi_e_face_op->SpMat().NumNonZeroElems();
				if (pml_upwind) {
					collectPMLUpwindDriverBlocks(
						blocks, psi_e_row, E, in_c, psi_zero_op->SpMat(),
						two_dir_pairs, psi_two_ptr_E, globalId, stats);
					if (isPMLOperatorAuditEnabled() && Mpi::WorldRank() == 0) {
						pmlAuditLog("PML_psiE_upwind_zero/two placed for in_c=" +
						            std::to_string(in_c));
					}
				}
				if (isPMLOperatorAuditEnabled() && Mpi::WorldRank() == 0) {
					pmlAuditLog("PML_psiE_driver row=" + std::to_string(psi_e_row) +
					            " Ey_col w=" + std::to_string(weight) +
					            " vol_nnz=" + std::to_string(psi_e_vol_op->SpMat().NumNonZeroElems()) +
					            " face_nnz=" + std::to_string(psi_e_face_op->SpMat().NumNonZeroElems()));
				}

				collectBlockPlacement(
					h_corr_op->SpMat(), blocks,
					globalId.offsets[H][h_comp]->rowStartOffset, psi_e_col, h_corr_sign);
				stats.correction_nnz += h_corr_op->SpMat().NumNonZeroElems();
			}

			for (Direction e_comp = X; e_comp <= Z; ++e_comp) {
				if (!pmlPsiHComponentActive(e_comp, stretch_dir)) {
					continue;
				}
				Direction in_c = X;
				double weight = 0.0;
				if (!pmlPsiHDriverCoupling(e_comp, stretch_dir, in_c, weight)) {
					continue;
				}

				const int psi_h_row = layout.psiHOffset(stretch_dir, e_comp);
				const int psi_h_col = 6 * blockSize + (psi_h_row - 6 * ndofs);

				collectBlockPlacement(
					psi_mass_op->SpMat(), blocks, psi_h_row, psi_h_col, psi_mass_sign);
				collectPMLComponentDriverBlocks(
					blocks, psi_h_row, H, in_c, stretch_dir, weight,
					psi_h_vol_op->SpMat(), psi_h_face_op->SpMat(),
					globalId, stats);
				stats.driver_volume_nnz += psi_h_vol_op->SpMat().NumNonZeroElems();
				stats.driver_face_nnz += psi_h_face_op->SpMat().NumNonZeroElems();
				if (pml_upwind) {
					collectPMLUpwindDriverBlocks(
						blocks, psi_h_row, H, in_c, psi_zero_op->SpMat(),
						two_dir_pairs, psi_two_ptr_H, globalId, stats);
				}

				collectBlockPlacement(
					e_corr_op->SpMat(), blocks,
					globalId.offsets[E][e_comp]->rowStartOffset, psi_h_col, e_corr_sign);
				stats.correction_nnz += e_corr_op->SpMat().NumNonZeroElems();
			}
		}
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildPMLOperator(
		const PMLAuxLayout& layout)
	{
		if (layout.nAux() == 0) {
			return nullptr;
		}

		const PMLProfileData* profiles = pd_.model.getPMLProfileData();
		if (!profiles) {
			throw std::runtime_error("buildPMLOperator requires initialized PML profile data.");
		}

		const int ndofs = fes_.GetNDofs();
		const int nbrDofs = getAdditionalDofs();
		const int blockSize = ndofs + nbrDofs;
		const int globalRows = 6 * ndofs + layout.nAux();
		const int globalCols = 6 * blockSize + layout.nAux();

		mfem::Array<int> pml_marker = pd_.model.buildPMLVolumeMarker();
		auto MInvScalar = buildPMLScalarInverseMassSubOperator<ParBilinearForm>(
			pml_marker, profiles->feOrder() + 1);
		auto MInvMaxwell = buildMaxwellInverseMassMatrixOperator<ParBilinearForm>();

		PMLOperatorAssemblyStats stats;
		std::vector<CSRBlockPlacement> blocks;
		collectPMLOperatorBlocks(
			blocks, layout, *profiles, pml_marker, MInvScalar, MInvMaxwell, stats);

		auto res = mergeBlocksToCSR(blocks, globalRows, globalCols);
		blocks.clear();
		res->Threshold(1e-8);

		if (Mpi::WorldRank() == 0) {
			std::cout << "[PML] PMLOperator_ assembled: " << globalRows << " x " << globalCols
			          << ", nnz=" << res->NumNonZeroElems()
			          << ", stretch_dirs=" << layout.numStretchDirections()
			          << ", driver_vol_nnz=" << stats.driver_volume_nnz
			          << ", driver_face_nnz=" << stats.driver_face_nnz
			          << ", driver_upwind_nnz=" << stats.driver_upwind_nnz
			          << ", upwind_alpha=" << pd_.opts.alpha
			          << ", corr_nnz=" << stats.correction_nnz
			          << ", sign_test=" << pmlSignTestModeName(getPMLSignTestMode())
			          << ", mult_sign=" << getPMLOperatorMultSign()
			          << std::endl;
		}

		if (pd_.opts.export_evolution_operator && Mpi::WorldSize() == 1) {
			std::filesystem::path export_dir =
				std::filesystem::path("Exports") / "Operators" / pd_.model.meshName_;
			if (!std::filesystem::exists(export_dir)) {
				std::filesystem::create_directories(export_dir);
			}
			std::filesystem::path file_path =
				export_dir / (pd_.model.meshName_ + "_pml.csr");
			std::ofstream ofs(file_path);
			if (ofs.is_open()) {
				res->PrintCSR2(ofs);
				ofs.close();
				if (Mpi::WorldRank() == 0) {
					std::cout << "PML operator exported to " << file_path << std::endl;
				}
			}
		}

		if (isPMLOperatorAuditEnabled()) {
			auditPMLOperatorCurl(layout, *res, pml_marker, *profiles);
		}

		return res;
	}

}