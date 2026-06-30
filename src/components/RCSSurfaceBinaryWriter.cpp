#include "RCSSurfaceBinaryWriter.h"

#include "components/FarField.h"

#include <filesystem>

namespace maxwell {

using namespace mfem;

RCSSurfaceBinaryWriter::RCSSurfaceBinaryWriter(
    const std::vector<int>& tags,
    const DG_FECollection* fec,
    ParFiniteElementSpace& parentFes,
    Fields<ParFiniteElementSpace, ParGridFunction>& globalFields,
    const std::string& outputRankPath)
    : submesher_(*parentFes.GetMesh(), parentFes, buildSurfaceMarker(tags, parentFes)),
      surfaceFes_(std::make_unique<FiniteElementSpace>(submesher_.getSubMesh(), fec)),
      surfaceFields_(*surfaceFes_),
      globalFields_(globalFields),
      transferMaps_(globalFields_, surfaceFields_)
{
    auto* mesh = submesher_.getSubMesh();
    spaceDim_ = mesh->SpaceDimension();
    numDofs_ = surfaceFields_.get(E, X).Size();

    std::filesystem::create_directories(outputRankPath);

    const auto elemOrder = parentFes.GetMesh()->GetElementTransformation(0)->Order();
    mesh->SetCurvature(elemOrder);
    mesh->Save(outputRankPath + "/mesh");

    dataFile_.open(outputRankPath + "/surface_data.bin", std::ios::binary);
    if (!dataFile_) {
        throw std::runtime_error("Could not open surface_data.bin for writing: " + outputRankPath);
    }

    const int numBdr = mesh->GetNBE();
    const auto bdrMarker = getNearToFarFieldMarker(mesh->bdr_attributes.Max());

    int totalQuadPts = 0;
    for (int be = 0; be < numBdr; ++be) {
        if (bdrMarker[mesh->GetBdrAttribute(be) - 1] != 1) continue;
        auto* Tr = mesh->GetBdrFaceTransformations(be);
        if (!Tr) continue;
        const auto& el = *surfaceFes_->GetFE(Tr->Elem1No);
        const auto* ir = &IntRules.Get(Tr->GetGeometryType(), 2 * (el.GetOrder() + Tr->Order()));
        totalQuadPts += ir->GetNPoints();
    }

    const int basisType = fec->GetBasisType();
    const int32_t header[5] = {
        static_cast<int32_t>(spaceDim_),
        static_cast<int32_t>(numDofs_),
        static_cast<int32_t>(numBdr),
        static_cast<int32_t>(totalQuadPts),
        static_cast<int32_t>(basisType)
    };
    dataFile_.write(reinterpret_cast<const char*>(header), sizeof(header));

    writeGeometry();
}

void RCSSurfaceBinaryWriter::writeGeometry()
{
    auto* mesh = submesher_.getSubMesh();
    const auto bdrMarker = getNearToFarFieldMarker(mesh->bdr_attributes.Max());

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
            const double face_weight = Tr->Weight();

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
}

void RCSSurfaceBinaryWriter::writeSnapshot(double time)
{
    transferMaps_.transferFields(globalFields_, surfaceFields_);

    dataFile_.write(reinterpret_cast<const char*>(&time), sizeof(double));

    for (auto ft : {E, H}) {
        for (auto d : {X, Y, Z}) {
            const auto& gf = surfaceFields_.get(ft, d);
            dataFile_.write(reinterpret_cast<const char*>(gf.GetData()),
                            numDofs_ * sizeof(double));
        }
    }
    dataFile_.flush();
}

} // namespace maxwell
