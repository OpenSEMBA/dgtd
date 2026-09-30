#include "SourcesManager.h"

#include "general/forall.hpp"

namespace maxwell {

namespace {

constexpr int kDirectSourceMetaStride = 2;
constexpr int kDirectSourceParamStride = 9;
constexpr double kPi = 3.14159265358979323846;

MFEM_HOST_DEVICE inline int metaIndex(int source_idx, int offset)
{
    return source_idx * kDirectSourceMetaStride + offset;
}

MFEM_HOST_DEVICE inline int paramIndex(int source_idx, int offset)
{
    return source_idx * kDirectSourceParamStride + offset;
}

MFEM_HOST_DEVICE inline double evalWaveform(int waveform,
                                            double spread,
                                            double mean,
                                            double frequency,
                                            double phase_delay,
                                            double time)
{
    const double arg = phase_delay - time - mean;
    const double envelope = exp(-(arg * arg) / (2.0 * spread * spread));
    if (waveform == static_cast<int>(DirectSourceWaveform::ModulatedGaussian)) {
        return envelope * cos(2.0 * kPi * frequency * arg);
    }
    return envelope;
}

} // namespace

void eval_tfsf_planewave_gpu(const mfem::Vector& dof_coords,
                             const mfem::Vector& tfsf_sign,
                             const mfem::Vector& source_params,
                             const mfem::Array<int>& source_meta,
                             bool apply_tfsf_sign,
                             double time,
                             int ndofs,
                             FieldGridFuncs& fields)
{
    for (int ft : {E, H}) {
        for (int d : {X, Y, Z}) {
            fields[ft][d] = 0.0;
        }
    }

    const int source_count = source_meta.Size() / kDirectSourceMetaStride;
    const double* coords = dof_coords.Read();
    const double* sign = tfsf_sign.Read();
    const double* params = source_params.Read();
    const int* meta = source_meta.Read();

    double* e0 = fields[E][X].Write();
    double* e1 = fields[E][Y].Write();
    double* e2 = fields[E][Z].Write();
    double* h0 = fields[H][X].Write();
    double* h1 = fields[H][Y].Write();
    double* h2 = fields[H][Z].Write();

    mfem::forall(ndofs, [=] MFEM_DEVICE(int i) {
        const double x = coords[3 * i + 0];
        const double y = coords[3 * i + 1];
        const double z = coords[3 * i + 2];
        const double scale = apply_tfsf_sign ? sign[i] : 1.0;

        double e_out[3] = {0.0, 0.0, 0.0};
        double h_out[3] = {0.0, 0.0, 0.0};

        for (int s = 0; s < source_count; ++s) {
            const int field_type = meta[metaIndex(s, 0)];
            const int waveform = meta[metaIndex(s, 1)];
            const double spread = params[paramIndex(s, 0)];
            const double mean = params[paramIndex(s, 1)];
            const double frequency = params[paramIndex(s, 2)];
            const double pol0 = params[paramIndex(s, 3)];
            const double pol1 = params[paramIndex(s, 4)];
            const double pol2 = params[paramIndex(s, 5)];
            const double dir0 = params[paramIndex(s, 6)];
            const double dir1 = params[paramIndex(s, 7)];
            const double dir2 = params[paramIndex(s, 8)];

            const double phase_delay =
                (x * dir0 + y * dir1 + z * dir2) / physicalConstants::speedOfLight;
            const double waveform_val =
                evalWaveform(waveform, spread, mean, frequency, phase_delay, time);

            double cross0;
            double cross1;
            double cross2;
            if (field_type == static_cast<int>(E)) {
                cross0 = dir1 * pol2 - dir2 * pol1;
                cross1 = dir2 * pol0 - dir0 * pol2;
                cross2 = dir0 * pol1 - dir1 * pol0;
                e_out[0] += waveform_val * pol0;
                e_out[1] += waveform_val * pol1;
                e_out[2] += waveform_val * pol2;
                h_out[0] += waveform_val * cross0;
                h_out[1] += waveform_val * cross1;
                h_out[2] += waveform_val * cross2;
            } else {
                cross0 = pol1 * dir2 - pol2 * dir1;
                cross1 = pol2 * dir0 - pol0 * dir2;
                cross2 = pol0 * dir1 - pol1 * dir0;
                h_out[0] += waveform_val * pol0;
                h_out[1] += waveform_val * pol1;
                h_out[2] += waveform_val * pol2;
                e_out[0] += waveform_val * cross0;
                e_out[1] += waveform_val * cross1;
                e_out[2] += waveform_val * cross2;
            }
        }

        e0[i] = scale * e_out[0];
        e1[i] = scale * e_out[1];
        e2[i] = scale * e_out[2];
        h0[i] = scale * h_out[0];
        h1[i] = scale * h_out[1];
        h2[i] = scale * h_out[2];
    });
}

} // namespace maxwell
