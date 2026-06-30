#include "PMLOperatorAudit.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace maxwell {

namespace {

bool envFlag(const char* name)
{
	const char* v = std::getenv(name);
	return v && v[0] == '1' && v[1] == '\0';
}

} // namespace

bool isPMLOperatorAuditEnabled()
{
	return envFlag("PML_OPERATOR_AUDIT");
}

bool shouldExportPMLOperatorAuditCSR()
{
	return isPMLOperatorAuditEnabled() || envFlag("PML_EXPORT_AUDIT_CSR");
}

bool shouldSkipGlobalOperatorInMult()
{
	return envFlag("PML_SKIP_GLOBAL");
}

bool shouldSkipPMLOperatorInMult()
{
	return envFlag("PML_SKIP_PML");
}

bool isPMLMultProbeEnabled()
{
	return envFlag("PML_MULT_PROBE");
}

void pmlAuditLog(const std::string& line)
{
	if (isPMLOperatorAuditEnabled()) {
		std::cout << "[PML audit] " << line << std::endl;
	}
}

double frobeniusBlockNorm(const mfem::SparseMatrix& mat,
                          int row0, int row1,
                          int col0, int col1)
{
	double sum = 0.0;
	for (int i = row0; i < row1; ++i) {
		const int n = mat.RowSize(i);
		const int* cols = mat.GetRowColumns(i);
		const double* vals = mat.GetRowEntries(i);
		for (int j = 0; j < n; ++j) {
			if (cols[j] >= col0 && cols[j] < col1) {
				sum += vals[j] * vals[j];
			}
		}
	}
	return std::sqrt(sum);
}

double frobeniusBlockNormOnRows(const mfem::SparseMatrix& mat,
                                const std::vector<int>& row_dofs,
                                int row_block_offset,
                                int col0, int col1)
{
	double sum = 0.0;
	for (int dof : row_dofs) {
		const int i = row_block_offset + dof;
		if (i < 0 || i >= mat.Height()) {
			continue;
		}
		const int n = mat.RowSize(i);
		const int* cols = mat.GetRowColumns(i);
		const double* vals = mat.GetRowEntries(i);
		for (int j = 0; j < n; ++j) {
			if (cols[j] >= col0 && cols[j] < col1) {
				sum += vals[j] * vals[j];
			}
		}
	}
	return std::sqrt(sum);
}

void exportAuditCSR(const mfem::SparseMatrix& mat,
                    const std::string& export_dir,
                    const std::string& filename)
{
	std::ofstream ofs(export_dir + "/" + filename);
	if (ofs.is_open()) {
		mat.PrintCSR2(ofs);
		ofs.close();
		pmlAuditLog("exported " + export_dir + "/" + filename +
		            " nnz=" + std::to_string(mat.NumNonZeroElems()));
	}
}

void printIntegratorBaseline1DTE()
{
	std::cout << "[PML audit] === 1D TE integrator baseline (active_axes X) ===\n"
	          << "[PML audit] Global curl chain Ey -> Hz:\n"
	          << "[PML audit]   collectGlobalDirectionalOperators: D_x on volume, "
	            "MInv[H]*(D_x) maps row H[z] <- col E[y] weight -1\n"
	          << "[PML audit]   collectGlobalOneNormalOperators: OneNormal({x}) on "
	            "altField(E), same H[z]<-E[y] SBP face flux weight -1\n"
	          << "[PML audit]   Zero/Two: coeff=upwind_alpha (0 => no contribution)\n"
	          << "[PML audit] PML psi chain (Gedney):\n"
	          << "[PML audit]   psi^E_{X,Z} dot = -alpha*psi + sigma*(D_x vol + SBP face) "
	            "on E[y]; pmlCurlDerivWeight(H,Z,E,Y,X)=+1\n"
	          << "[PML audit]   H[z] dot -= psi^E/kappa (PMLOperator correction block)\n"
	          << "[PML audit] Duplicate-curl risk: global applies full D_x in PML; "
	            "PMLOperator adds parallel sigma*D_x into psi driver.\n"
	          << "[PML audit] ================================================\n";
}

} // namespace maxwell
