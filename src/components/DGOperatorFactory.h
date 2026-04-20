#pragma once

#include "ProblemDescription.h"

#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <type_traits>

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

	/// Merge collected CSR block placements into a single finalized CSR SparseMatrix.
	/// Uses a two-pass marker technique (like MFEM's Add) to avoid LIL overhead.
	inline std::unique_ptr<SparseMatrix> mergeBlocksToCSR(
		std::vector<CSRBlockPlacement>& blocks,
		int globalRows, int globalCols)
	{
		// Pass 1: Count unique column entries per global row.
		std::vector<int> marker(globalCols, -1);
		int* C_i = mfem::Memory<int>(globalRows + 1);
		C_i[0] = 0;

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
			C_i[row + 1] = C_i[row] + nnz;
		}

		const int totalNNZ = C_i[globalRows];
		int* C_j = mfem::Memory<int>(totalNNZ);
		real_t* C_data = mfem::Memory<real_t>(totalNNZ);

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

		return std::make_unique<SparseMatrix>(C_i, C_j, C_data, globalRows, globalCols);
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

	private:
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
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
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
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
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
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
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
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
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
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
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
		GlobalIndices globalId(fes_.GetNDofs(), getAdditionalDofs(), true);
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
		return buildSourceFaceOperator(BdrCond::SGBC);
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
			this->template addGlobalSourceFaceIBFIOneNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
			this->template addGlobalSourceFaceIBFIZeroNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
			this->template addGlobalSourceFaceIBFITwoNormalOperators<ParBilinearForm>(res.get(), marker, MInv);
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

		auto threshold = 1e-20;
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

			std::filesystem::path file_path = export_dir / "SpatialEvolutionOperator.csr";

			std::ofstream ofs(file_path);
			if (!ofs.is_open())
			{
				throw std::runtime_error("Could not open file for writing: " + file_path.string());
			}

			res->PrintCSR2(ofs);
			ofs.close();

			std::cout << "Operator exported to " << file_path << std::endl;
		}

		return res;
	}

}