#include "wires/ThinWire.h"

#include "general/forall.hpp"

#include <algorithm>
#include <cmath>

namespace maxwell {
namespace {

MFEM_HOST_DEVICE
double lineTrace(const double* dofs, const double* shape, int nloc)
{
	double value = 0.0;
	for (int i = 0; i < nloc; ++i) {
		value += dofs[i] * shape[i];
	}
	return value;
}

}  // namespace

void ThinWireCoupler::advanceOnDevice(
	double tau,
	const mfem::Vector& eh,
	mfem::Vector* out,
	int ndofs,
	mfem::Vector& line_current,
	mfem::Vector& line_charge,
	bool copy_from_committed,
	double sample_time) const
{
	const int nline = line_current.Size();
	const int npieces = static_cast<int>(pieces_.size());
	if (nline == 0 || npieces == 0) {
		return;
	}
	const int nloc = order_ + 1;
	const int nvol = nvol_;
	const int es_stride = nloc * nvol;
	const int j_stride = nvol * nloc;
	const int mat_stride = nloc * nloc;

	if (copy_from_committed) {
		const double* src_current = current_.Read();
		const double* src_charge = charge_.Read();
		double* dst_current = line_current.Write();
		double* dst_charge = line_charge.Write();
		mfem::forall(nline, [=] MFEM_HOST_DEVICE (int i) {
			dst_current[i] = src_current[i];
			dst_charge[i] = src_charge[i];
		});
	}

	if (tau > 0.0) {
		const double* field = eh.Read();
		const int* vol_dofs = vol_dofs_.Read();
		const int* block_piece = block_piece_.Read();
		const int* offset = offset_.Read();
		const double* es = es_.Read();
		const double* mass_inv = mass_inv_.Read();
		double* drive = drive_.Write();
		const int nblocks = nblocks_;
		const int feed_count = feed_count_;
		const int feed_piece0 = feed_piece_[0];
		const int feed_piece1 = feed_piece_[1];
		const double feed_field0 = feedElectric(sample_time) * feed_per_length_[0];
		const double feed_field1 = feedElectric(sample_time) * feed_per_length_[1];
		mfem::forall(nline, [=] MFEM_HOST_DEVICE (int i) {
			drive[i] = 0.0;
		});
		if (nblocks > 0 && nvol > 0) {
			double* drive_acc = drive_.ReadWrite();
			mfem::forall(nblocks, [=] MFEM_HOST_DEVICE (int b) {
				const int p = block_piece[b];
				double load[16];
				for (int row = 0; row < nloc; ++row) {
					double sum = 0.0;
					for (int d = 0; d < 3; ++d) {
						const double* block = es + (b * 3 + d) * es_stride + row * nvol;
						for (int col = 0; col < nvol; ++col) {
							const int dof = vol_dofs[b * nvol + col];
							sum += block[col] * field[d * ndofs + dof];
						}
					}
					load[row] = sum;
				}
				const double* minv = mass_inv + p * mat_stride;
				const int piece_offset = offset[p];
				for (int row = 0; row < nloc; ++row) {
					double nodal = 0.0;
					for (int col = 0; col < nloc; ++col) {
						nodal += minv[row * nloc + col] * load[col];
					}
					::AtomicAdd(drive_acc[piece_offset + row], nodal);
				}
			});
		}
		double* drive_feed = drive_.ReadWrite();
		mfem::forall(npieces, [=] MFEM_HOST_DEVICE (int p) {
			if (feed_count > 0 && p == feed_piece0) {
				const int piece_offset = offset[p];
				for (int row = 0; row < nloc; ++row) {
					drive_feed[piece_offset + row] += feed_field0;
				}
			}
			if (feed_count > 1 && p == feed_piece1) {
				const int piece_offset = offset[p];
				for (int row = 0; row < nloc; ++row) {
					drive_feed[piece_offset + row] += feed_field1;
				}
			}
		});

		const int steps = std::max(1, static_cast<int>(std::ceil(tau / substep_)));
		const double step = tau / static_cast<double>(steps);
		const double* stiffness = stiffness_.Read();
		const double* shape_left = shape_left_.Read();
		const double* shape_right = shape_right_.Read();
		const double* v2 = v2_.Read();
		const double* inductance = inductance_.Read();
		const double* resistance = resistance_.Read();
		const int* left_neighbor = left_neighbor_.Read();
		const int* right_neighbor = right_neighbor_.Read();
		const int* left_pec = left_pec_.Read();
		const int* right_pec = right_pec_.Read();
		const double* drive_read = drive_.Read();

		for (int sub = 0; sub < steps; ++sub) {
			const double* current = line_current.Read();
			const double* charge = line_charge.Read();
			double* d_current = d_current_.Write();
			double* d_charge = d_charge_.Write();
			mfem::forall(npieces, [=] MFEM_HOST_DEVICE (int p) {
				const int piece_offset = offset[p];
				double line_i[16];
				double line_q[16];
				for (int i = 0; i < nloc; ++i) {
					line_i[i] = current[piece_offset + i];
					line_q[i] = charge[piece_offset + i];
				}
				const double* stiff = stiffness + p * mat_stride;
				const double* minv = mass_inv + p * mat_stride;
				double stiff_i[16];
				double stiff_q[16];
				double vol_i[16];
				double vol_q[16];
				for (int row = 0; row < nloc; ++row) {
					double sum_i = 0.0;
					double sum_q = 0.0;
					for (int col = 0; col < nloc; ++col) {
						const double weight = stiff[row * nloc + col];
						sum_i += weight * line_i[col];
						sum_q += weight * line_q[col];
					}
					stiff_i[row] = sum_i;
					stiff_q[row] = sum_q;
				}
				for (int row = 0; row < nloc; ++row) {
					double sum_i = 0.0;
					double sum_q = 0.0;
					for (int col = 0; col < nloc; ++col) {
						const double weight = minv[row * nloc + col];
						sum_i += weight * stiff_i[col];
						sum_q += weight * stiff_q[col];
					}
					vol_i[row] = sum_i;
					vol_q[row] = sum_q;
				}

				const double piece_v2 = v2[p];
				const double piece_inductance = inductance[p];
				const double piece_resistance = resistance[p];
				const double speed = sqrt(piece_v2);
				double d_i[16];
				double d_q[16];
				for (int i = 0; i < nloc; ++i) {
					d_i[i] = -piece_v2 * vol_q[i] - (piece_resistance / piece_inductance) * line_i[i];
					d_q[i] = -vol_i[i];
					d_i[i] += drive_read[piece_offset + i] / piece_inductance;
				}

				for (int end = 0; end < 2; ++end) {
					const bool right = end == 1;
					const double* shape = right ? shape_right : shape_left;
					const double own_i = lineTrace(line_i, shape, nloc);
					const double own_q = lineTrace(line_q, shape, nloc);
					const int neighbor = right ? right_neighbor[p] : left_neighbor[p];
					const int pec = right ? right_pec[p] : left_pec[p];
					const double alpha_out = right
						? 0.5 * (own_q + own_i / speed)
						: 0.5 * (own_q - own_i / speed);
					double alpha_in = 0.0;
					if (neighbor >= 0) {
						const int neighbor_offset = offset[neighbor];
						double neighbor_i[16];
						double neighbor_q[16];
						for (int i = 0; i < nloc; ++i) {
							neighbor_i[i] = current[neighbor_offset + i];
							neighbor_q[i] = charge[neighbor_offset + i];
						}
						const double* neighbor_shape = right ? shape_left : shape_right;
						const double exterior_i = lineTrace(neighbor_i, neighbor_shape, nloc);
						const double exterior_q = lineTrace(neighbor_q, neighbor_shape, nloc);
						alpha_in = right
							? 0.5 * (exterior_q - exterior_i / speed)
							: 0.5 * (exterior_q + exterior_i / speed);
					} else if (pec != 0) {
						alpha_in = -alpha_out;
					} else {
						alpha_in = alpha_out;
					}
					double hat_i = 0.0;
					double hat_q = 0.0;
					if (right) {
						hat_i = speed * (alpha_out - alpha_in);
						hat_q = alpha_out + alpha_in;
					} else {
						hat_i = speed * (alpha_in - alpha_out);
						hat_q = alpha_in + alpha_out;
					}
					const double normal = right ? 1.0 : -1.0;
					const double flux_i = normal * piece_v2 * (own_q - hat_q);
					const double flux_q = normal * (own_i - hat_i);
					double lift_i[16];
					double lift_q[16];
					for (int i = 0; i < nloc; ++i) {
						lift_i[i] = shape[i] * flux_i;
						lift_q[i] = shape[i] * flux_q;
					}
					for (int row = 0; row < nloc; ++row) {
						double add_i = 0.0;
						double add_q = 0.0;
						for (int col = 0; col < nloc; ++col) {
							const double weight = minv[row * nloc + col];
							add_i += weight * lift_i[col];
							add_q += weight * lift_q[col];
						}
						d_i[row] += add_i;
						d_q[row] += add_q;
					}
				}

				for (int i = 0; i < nloc; ++i) {
					d_current[piece_offset + i] = d_i[i];
					d_charge[piece_offset + i] = d_q[i];
				}
			});

			const double* rate_current = d_current_.Read();
			const double* rate_charge = d_charge_.Read();
			double* next_current = line_current.ReadWrite();
			double* next_charge = line_charge.ReadWrite();
			mfem::forall(nline, [=] MFEM_HOST_DEVICE (int i) {
				next_current[i] += step * rate_current[i];
				next_charge[i] += step * rate_charge[i];
			});
		}
	}

	if (out != nullptr && nblocks_ > 0 && nvol > 0) {
		const double* current = line_current.Read();
		const int* vol_dofs = vol_dofs_.Read();
		const int* block_piece = block_piece_.Read();
		const int* offset = offset_.Read();
		const double* j_mat = j_mat_.Read();
		double* residual = out->ReadWrite();
		const double return_scale = electric_return_scale_;
		const int nblocks = nblocks_;
		mfem::forall(nblocks, [=] MFEM_HOST_DEVICE (int b) {
			const int p = block_piece[b];
			const int piece_offset = offset[p];
			for (int d = 0; d < 3; ++d) {
				const double* block = j_mat + (b * 3 + d) * j_stride;
				for (int row = 0; row < nvol; ++row) {
					double sum = 0.0;
					for (int col = 0; col < nloc; ++col) {
						sum += block[row * nloc + col] * current[piece_offset + col];
					}
					const int dof = vol_dofs[b * nvol + row];
					::AtomicAdd(residual[d * ndofs + dof], return_scale * sum);
				}
			}
		});
	}
}

}  // namespace maxwell
