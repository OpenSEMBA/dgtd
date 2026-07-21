#include "RCSSurfaceExporter.h"
#include "components/FarField.h"
#include "math/PhysicalConstants.h"
#include "solver/ProbesManager.h"

namespace maxwell {

using namespace mfem;

RCSSurfaceExporter::RCSSurfaceExporter(
    const RCSSurfaceProbe& probe,
    const DG_FECollection* fec,
    ParFiniteElementSpace& parentFes,
    Fields<ParFiniteElementSpace, ParGridFunction>& globalFields,
    const std::string& caseName)
    : submesher_(*parentFes.GetMesh(), parentFes, buildSurfaceMarker(probe.tags, parentFes)),
      globalFields_(globalFields),
      expSteps_(probe.expSteps)
{
    if (!submesher_.hasLocalSurface()) {
        spaceDim_ = parentFes.GetMesh()->SpaceDimension();
        return;
    }

    hasLocalSurface_ = true;
    surfaceFes_ = std::make_unique<FiniteElementSpace>(submesher_.getSubMesh(), fec);
    surfaceFields_ = std::make_unique<Fields<FiniteElementSpace, GridFunction>>(*surfaceFes_);
    transferMaps_ = std::make_unique<TransferMaps>(globalFields, *surfaceFields_);
    initParentToSurfaceMap();

    auto* mesh = submesher_.getSubMesh();
    spaceDim_ = mesh->SpaceDimension();
    numDofs_ = surfaceFields_->get(E, X).Size();

    std::string base = "Exports/" + getRunModeTag() + "/" + caseName + "/RCSSurface/" + probe.name;
    std::filesystem::create_directories(base);
    outputPath_ = base + "/rank" + std::to_string(Mpi::WorldRank());
    std::filesystem::create_directories(outputPath_);

    auto elemOrder = parentFes.GetMesh()->GetElementTransformation(0)->Order();
    mesh->SetCurvature(elemOrder);
    mesh->Save(outputPath_ + "/mesh");

    dataFile_.open(outputPath_ + "/surface_data.bin", std::ios::binary);

    // Count total quadrature points on NTF boundary faces.
    int numBdr = mesh->GetNBE();
    auto bdrMarker = getNearToFarFieldMarker(mesh->bdr_attributes.Max());

    int totalQuadPts = 0;
    for (int be = 0; be < numBdr; ++be) {
        if (bdrMarker[mesh->GetBdrAttribute(be) - 1] != 1) continue;
        auto* Tr = mesh->GetBdrFaceTransformations(be);
        if (!Tr) continue;
        const auto& el = *surfaceFes_->GetFE(Tr->Elem1No);
        const auto* ir = &IntRules.Get(Tr->GetGeometryType(), 2 * (el.GetOrder() + Tr->Order()));
        totalQuadPts += ir->GetNPoints();
    }

    int basisType = fec->GetBasisType();
    int32_t header[5] = {
        static_cast<int32_t>(spaceDim_),
        static_cast<int32_t>(numDofs_),
        static_cast<int32_t>(numBdr),
        static_cast<int32_t>(totalQuadPts),
        static_cast<int32_t>(basisType)
    };
    dataFile_.write(reinterpret_cast<const char*>(header), sizeof(header));

    writeGeometry();
    transferFields();
}

void RCSSurfaceExporter::initParentToSurfaceMap()
{
    auto* surface_sm = static_cast<SubMesh*>(surfaceFes_->GetMesh());
    Array<int> raw_map;
    SubMeshUtils::BuildVdofToVdofMap(*surfaceFes_,
                                     *globalFields_.get(E, X).ParFESpace(),
                                     surface_sm->GetFrom(),
                                     surface_sm->GetParentElementIDMap(),
                                     raw_map);

    parent_dof_ids_.SetSize(raw_map.Size());
    parent_dof_signs_.SetSize(raw_map.Size());
    for (int i = 0; i < raw_map.Size(); ++i) {
        real_t sign = 1.0;
        parent_dof_ids_[i] = FiniteElementSpace::DecodeDof(raw_map[i], sign);
        parent_dof_signs_[i] = sign;
    }
}

void RCSSurfaceExporter::transferFields()
{
    if (!hasLocalSurface_) return;

#ifdef SEMBA_DGTD_ENABLE_CUDA
    if (mfem::Device::Allows(mfem::Backend::CUDA)) {
        rcs_gather_surface_fields_gpu(parent_dof_ids_,
                                      parent_dof_signs_,
                                      globalFields_,
                                      *surfaceFields_,
                                      numDofs_);
        return;
    }
#endif
    transferMaps_->transferFields(globalFields_, *surfaceFields_);
}

void RCSSurfaceExporter::writeGeometry()
{
    if (!hasLocalSurface_) return;

    auto* mesh = submesher_.getSubMesh();
    auto bdrMarker = getNearToFarFieldMarker(mesh->bdr_attributes.Max());

    std::vector<double> positions;
    std::vector<double> normals;
    std::vector<double> weights;

    for (int be = 0; be < mesh->GetNBE(); ++be) {
        if (bdrMarker[mesh->GetBdrAttribute(be) - 1] != 1) continue;

        auto* Tr = mesh->GetBdrFaceTransformations(be);
        if (!Tr) continue;

        const auto& el = *surfaceFes_->GetFE(Tr->Elem1No);
        const auto* ir = &IntRules.Get(Tr->GetGeometryType(), 2 * (el.GetOrder() + Tr->Order()));

        for (int q = 0; q < ir->GetNPoints(); ++q) {
            const auto& ip = ir->IntPoint(q);
            Tr->SetAllIntPoints(&ip);

            Vector phys_pt;
            Tr->Face->Transform(ip, phys_pt);
            for (int d = 0; d < spaceDim_; ++d) {
                positions.push_back(phys_pt(d));
            }

            Vector ortho(el.GetDim());
            CalcOrtho(Tr->Jacobian(), ortho);
            double face_weight = Tr->Weight();

            Vector outward(3);
            outward = 0.0;
            for (int d = 0; d < ortho.Size(); ++d) {
                outward[d] = ortho[d] / face_weight;
            }
            for (int d = 0; d < 3; ++d) {
                normals.push_back(outward[d]);
            }

            weights.push_back(ip.weight * face_weight);
        }
    }

    dataFile_.write(reinterpret_cast<const char*>(positions.data()),
                    positions.size() * sizeof(double));
    dataFile_.write(reinterpret_cast<const char*>(normals.data()),
                    normals.size() * sizeof(double));
    dataFile_.write(reinterpret_cast<const char*>(weights.data()),
                    weights.size() * sizeof(double));
    dataFile_.flush();
    geometryWritten_ = true;
}

void RCSSurfaceExporter::writeSnapshot(double time)
{
    if (!hasLocalSurface_) return;

    double t = time;
    dataFile_.write(reinterpret_cast<const char*>(&t), sizeof(double));

    for (auto ft : {E, H}) {
        for (auto d : {X, Y, Z}) {
            auto& gf = surfaceFields_->get(ft, d);
#ifdef SEMBA_DGTD_ENABLE_CUDA
            if (mfem::Device::Allows(mfem::Backend::CUDA)) {
                gf.HostRead();
            }
#endif
            dataFile_.write(reinterpret_cast<const char*>(gf.GetData()),
                            numDofs_ * sizeof(double));
        }
    }
    dataFile_.flush();
}

void RCSSurfaceExporter::write(double time, int cycle, double finalTime)
{
    if (!hasLocalSurface_) return;

    bool atEnd = std::abs(time - finalTime) < 1e-8;
    if (!atEnd && cycle % expSteps_ != 0) {
        return;
    }

    transferFields();
    writeSnapshot(time);
}

} // namespace maxwell
