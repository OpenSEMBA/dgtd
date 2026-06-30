#pragma once

#include "Types.h"

#include <mfem.hpp>

#include <string>
#include <vector>

namespace maxwell {

class PMLAuxLayout;

bool isPMLOperatorAuditEnabled();
bool shouldExportPMLOperatorAuditCSR();
bool shouldSkipGlobalOperatorInMult();
bool shouldSkipPMLOperatorInMult();
bool isPMLMultProbeEnabled();

void pmlAuditLog(const std::string& line);

/// Frobenius norm of entries with row in [row0,row1) and col in [col0,col1).
double frobeniusBlockNorm(const mfem::SparseMatrix& mat,
                          int row0, int row1,
                          int col0, int col1);

/// Frobenius norm restricted to rows in dof_list (cols in [col0,col1)).
double frobeniusBlockNormOnRows(const mfem::SparseMatrix& mat,
                                const std::vector<int>& row_dofs,
                                int row_block_offset,
                                int col0, int col1);

void exportAuditCSR(const mfem::SparseMatrix& mat,
                    const std::string& export_dir,
                    const std::string& filename);

void printIntegratorBaseline1DTE();

} // namespace maxwell
