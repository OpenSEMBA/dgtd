#pragma once

#include <fstream>
#include <string>
#include <filesystem>

#include "components/Probes.h"
#include "components/SubMesher.h"

namespace maxwell {

/**
 * Exports surface field data needed for RCS post-processing in a compact
 * binary format.  One file per MPI rank is produced containing:
 *
 *   HEADER
 *   ------
 *   int32  spaceDimension
 *   int32  numDofs            (number of DOFs per scalar field on the surface)
 *   int32  numBdrElements     (boundary elements on the NTF surface)
 *   int32  numQuadPoints      (total quadrature points across all bdr elements)
 *
 *   GEOMETRY (written once)
 *   -----------------------
 *   double[numQuadPoints * spaceDimension]   quadrature point positions
 *   double[numQuadPoints * 3]                outward unit normals (always 3-comp)
 *   double[numQuadPoints]                    quadrature weights (include Jacobian)
 *
 *   FIELD SNAPSHOTS (appended each export step)
 *   --------------------------------------------
 *   double  time
 *   double[numDofs]  Ex
 *   double[numDofs]  Ey
 *   double[numDofs]  Ez
 *   double[numDofs]  Hx
 *   double[numDofs]  Hy
 *   double[numDofs]  Hz
 *
 * The companion mesh is also saved once for reference / verification.
 */
class RCSSurfaceExporter {
public:
    RCSSurfaceExporter(
        const RCSSurfaceProbe& probe,
        const mfem::DG_FECollection* fec,
        mfem::ParFiniteElementSpace& parentFes,
        Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& globalFields,
        const std::string& caseName);

    void write(double time, int cycle, double finalTime);

    const std::string& getOutputPath() const { return outputPath_; }

private:
    void transferFields();
    void writeGeometry();
    void writeSnapshot(double time);

    NearToFarFieldSubMesher submesher_;
    std::unique_ptr<mfem::FiniteElementSpace> surfaceFes_;
    Fields<mfem::FiniteElementSpace, mfem::GridFunction> surfaceFields_;
    Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& globalFields_;

    mfem::TransferMap tMapEx_, tMapEy_, tMapEz_;
    mfem::TransferMap tMapHx_, tMapHy_, tMapHz_;

    std::string outputPath_;
    std::ofstream dataFile_;
    int numDofs_;
    int spaceDim_;
    int expSteps_;
    bool geometryWritten_{false};
};

} // namespace maxwell
