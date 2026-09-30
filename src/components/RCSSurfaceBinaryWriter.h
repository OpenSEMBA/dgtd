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
