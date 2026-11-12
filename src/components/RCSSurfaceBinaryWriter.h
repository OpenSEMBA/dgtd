#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "components/SubMesher.h"
#include "evolution/Fields.h"

namespace maxwell {

/**
 * Writes rcssurface probe binary export (rank0/mesh + surface_data.bin)
 * from full-order field snapshots.  Format matches RCSSurfaceExporter.
 */
class RCSSurfaceBinaryWriter {
public:
    RCSSurfaceBinaryWriter(
        const std::vector<int>& tags,
        const mfem::DG_FECollection* fec,
        mfem::ParFiniteElementSpace& parentFes,
        Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& globalFields,
        const std::string& outputRankPath);

    void writeSnapshot(double time);

    /// Scalar surface DOFs per field component (Ex, Ey, Ez, Hx, Hy, Hz).
    int numSurfaceDofs() const { return numDofs_; }

    /// Parent volume DOF index and sign for each surface DOF of one component.
    /// The same map applies to every component, with offset c * ndofs.
    void surfaceParentMap(mfem::Array<int>& ids, std::vector<double>& signs) const;

    /// Write one snapshot already ordered as Ex, Ey, Ez, Hx, Hy, Hz.
    /// `packed` has 6 * numSurfaceDofs() doubles. Does not flush.
    void writePackedSnapshot(double time, const double* packed);

    /// Pack Ex, Ey, Ez, Hx, Hy, Hz after the same transfer writeSnapshot uses.
    void copyTransferredSurface(double* packed);

    void flush();

private:
    void writeGeometry();

    NearToFarFieldSubMesher submesher_;
    std::unique_ptr<mfem::FiniteElementSpace> surfaceFes_;
    Fields<mfem::FiniteElementSpace, mfem::GridFunction> surfaceFields_;
    Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& globalFields_;
    TransferMaps transferMaps_;

    std::ofstream dataFile_;
    int numDofs_ = 0;
    int spaceDim_ = 0;
};

} // namespace maxwell
