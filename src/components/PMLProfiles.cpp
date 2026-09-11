#include "PMLProfiles.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace maxwell {

namespace {

double designSigmaMax(double L, const PMLProperties& props)
{
	if (L <= 0.0 || props.target_reflection <= 0.0 || props.target_reflection >= 1.0) {
		return 0.0;
	}
	const int m = props.grading_order;
	return -(static_cast<double>(m) + 1.0) * std::log(props.target_reflection) / (2.0 * L);
}

/// ∫_{r_inner}^{r_inner+ρ} σ_r(r') dr' for power-law / constant grading.
double integratedRadialSigma(double rho, double L, double sigma_max, int m)
{
	if (L <= 0.0 || rho <= 0.0 || sigma_max <= 0.0) {
		return 0.0;
	}
	const double xi = std::clamp(rho / L, 0.0, 1.0);
	// σ_r = σ_max (ρ'/L)^m  ⇒  Σ = σ_max L/(m+1) ξ^{m+1}  (m=0: Σ = σ_max ρ).
	return sigma_max * L / (static_cast<double>(m) + 1.0) *
	       std::pow(xi, static_cast<double>(m) + 1.0);
}

/// Cylindrical SC metric: s_θ = r̃/r ⇒ σ_θ = Σ(r)/r with Σ = ∫ σ_r dr'.
double cylindricalSigmaTheta(double r, double rho, double L, double sigma_max, int m)
{
	constexpr double r_eps = 1e-14;
	if (r < r_eps) {
		return 0.0;
	}
	return integratedRadialSigma(rho, L, sigma_max, m) / r;
}

void evaluateStretchProfiles(
	double rho, double L, const PMLProperties& props, PMLDirectionProfiles& out)
{
	out.depth = rho;
	out.sigma = 0.0;
	out.kappa = 1.0;
	out.alpha = 0.0;
	if (L <= 0.0) {
		return;
	}

	const double xi = std::clamp(rho / L, 0.0, 1.0);
	const int m = props.grading_order;
	const double sigma_max = designSigmaMax(L, props);

	if (m == 0) {
		out.sigma = sigma_max;
		out.kappa = props.kappa_max;
		out.alpha = props.alpha_max;
	} else {
		const double xi_m = std::pow(xi, static_cast<double>(m));
		out.sigma = sigma_max * xi_m;
		out.kappa = 1.0 + (props.kappa_max - 1.0) * xi_m;
		out.alpha = props.alpha_max * xi_m;
	}
}

int dominantAxis(const mfem::Vector& normal, int mesh_dim)
{
	int best = 0;
	double best_val = 0.0;
	for (int d = 0; d < mesh_dim; ++d) {
		const double val = std::abs(normal(d));
		if (val > best_val) {
			best_val = val;
			best = d;
		}
	}
	return best;
}

double distanceFromCenter(const mfem::Vector& x, const std::array<double, 3>& c, int dim)
{
	double r2 = 0.0;
	for (int d = 0; d < dim; ++d) {
		const double dx = x(d) - c[static_cast<size_t>(d)];
		r2 += dx * dx;
	}
	return std::sqrt(r2);
}

} // namespace

PMLProfileData::PMLProfileData(
	mfem::Mesh& mesh, const std::vector<PMLProperties>& regions, int fe_order)
	: fe_order_(fe_order), mesh_dim_(mesh.Dimension())
{
	if (regions.empty()) {
		return;
	}
	regions_ = regions;

	buildAttributeMaps(mesh, regions);
	buildInterfaceData(mesh, regions);
	buildRadialRegionData(mesh, regions);
	buildElementProfiles(mesh, regions, fe_order);
}

const PMLElementProfiles* PMLProfileData::getElementProfiles(int el) const
{
	if (el < 0 || el >= static_cast<int>(el_to_profile_index_.size())) {
		return nullptr;
	}
	const int index = el_to_profile_index_[el];
	if (index < 0) {
		return nullptr;
	}
	return &element_profiles_[index];
}

const PMLDirectionProfiles* PMLProfileData::getDirectionProfileAtIP(
	int el, const mfem::IntegrationPoint& ip, Direction stretch_dir) const
{
	const PMLElementProfiles* ep = getElementProfiles(el);
	if (!ep || stretch_dir < 0 || stretch_dir >= 3) {
		return nullptr;
	}
	if (stretch_dir >= 3 || ep->region_index < 0 ||
	    ep->region_index >= static_cast<int>(regions_.size())) {
		return nullptr;
	}

	const PMLProperties& props = regions_[ep->region_index];
	if (props.active_axes.count(stretch_dir) == 0) {
		return nullptr;
	}

	const int nqp = static_cast<int>(ep->qp_profiles.size());
	if (nqp == 0) {
		return nullptr;
	}

	const int order = fe_order_ + 1;
	const mfem::IntegrationRule& ir =
		mfem::IntRules.Get(mfem::Geometry::SEGMENT, order);
	const mfem::IntegrationRule* rule = &ir;
	if (static_cast<int>(ir.GetNPoints()) != nqp) {
		const mfem::IntegrationRule& ir2 =
			mfem::IntRules.Get(mfem::Geometry::TRIANGLE, order);
		if (static_cast<int>(ir2.GetNPoints()) == nqp) {
			rule = &ir2;
		}
	}

	int iq = -1;
	for (int i = 0; i < rule->GetNPoints() && i < nqp; ++i) {
		const mfem::IntegrationPoint& q = (*rule)[i];
		if (std::abs(q.x - ip.x) < 1e-12 && std::abs(q.y - ip.y) < 1e-12 &&
		    std::abs(q.z - ip.z) < 1e-12) {
			iq = i;
			break;
		}
	}
	if (iq < 0 && nqp == 1) {
		iq = 0;
	}
	if (iq < 0 || iq >= nqp) {
		return nullptr;
	}
	return &ep->qp_profiles[iq][stretch_dir];
}

void PMLProfileData::evaluateAtTransform(
	mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip,
	Direction stretch_dir, PMLDirectionProfiles& out) const
{
	out = PMLDirectionProfiles{};
	// Attribute-based region lookup — required for MPI. Profiles may be built on
	// the full serial mesh (global interface / L), while ParBilinearForm passes
	// ElementTransformations whose ElementNo is the *local* ParMesh index.
	const int attr = T.Attribute;
	if (attr <= 0 || attr > static_cast<int>(is_pml_attr_.size()) ||
	    !is_pml_attr_[attr - 1] || stretch_dir < 0 || stretch_dir >= 3) {
		return;
	}

	const int region_index = attr_region_index_[attr - 1];
	if (region_index < 0 || region_index >= static_cast<int>(regions_.size())) {
		return;
	}

	const PMLProperties& props = regions_[region_index];
	if (props.active_axes.count(stretch_dir) == 0) {
		return;
	}

	const int dim = T.GetSpaceDim();
	T.SetIntPoint(&ip);
	mfem::Vector x(dim);
	T.Transform(ip, x);

	double rho = 0.0;
	double L = 0.0;
	if (props.stretch_mode == PMLStretchMode::Radial &&
	    region_index < static_cast<int>(radial_.size()) &&
	    radial_[region_index].active) {
		rho = depthRadial(x, region_index);
		L = radial_[region_index].thickness();
	} else {
		rho = depthAlongAxis(x, stretch_dir);
		L = region_max_depth_[region_index][stretch_dir];
	}
	evaluateStretchProfiles(rho, L, props, out);
}

bool PMLProfileData::hasUniaxialRadialRegion() const
{
	for (const auto& props : regions_) {
		if (props.uniaxial_radial) {
			return true;
		}
	}
	return false;
}

bool PMLProfileData::evaluateRotatedTensorAtTransform(
	mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip,
	int tensor_kind, double T_xyz[3][3]) const
{
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			T_xyz[i][j] = 0.0;
		}
	}

	const int attr = T.Attribute;
	if (attr <= 0 || attr > static_cast<int>(is_pml_attr_.size()) ||
	    !is_pml_attr_[attr - 1]) {
		return false;
	}
	const int region_index = attr_region_index_[attr - 1];
	if (region_index < 0 || region_index >= static_cast<int>(regions_.size())) {
		return false;
	}
	const PMLProperties& props = regions_[region_index];
	if (!props.uniaxial_radial) {
		return false;
	}
	if (region_index >= static_cast<int>(radial_.size()) || !radial_[region_index].active) {
		return false;
	}

	const int dim = T.GetSpaceDim();
	T.SetIntPoint(&ip);
	mfem::Vector x(dim);
	T.Transform(ip, x);

	const RadialRegionData& rd = radial_[region_index];
	double dx[3] = {0.0, 0.0, 0.0};
	double r2 = 0.0;
	for (int d = 0; d < dim; ++d) {
		dx[d] = x(d) - rd.center[static_cast<size_t>(d)];
		r2 += dx[d] * dx[d];
	}
	const double r = std::sqrt(r2);
	constexpr double r_eps = 1e-14;
	if (r < r_eps) {
		return false;
	}

	const double rho = std::max(0.0, r - rd.r_inner);
	const double L = rd.thickness();
	PMLDirectionProfiles stretch;
	evaluateStretchProfiles(rho, L, props, stretch);
	const double sig_r = stretch.sigma;
	const double sigma_max = designSigmaMax(L, props);
	const double sig_theta =
		cylindricalSigmaTheta(r, rho, L, sigma_max, props.grading_order);
	// κ≡1 MVP for uniaxial_radial. Cylindrical SC: σ = (σ_r, σ_θ, 0) with
	// σ_θ = Σ/r, Σ = ∫_{r_inner}^r σ_r dr' (not locally uniaxial σ_θ=0).
	const double kap[3] = {1.0, 1.0, 1.0};
	const double sig[3] = {sig_r, sig_theta, 0.0};

	double prin[3];
	for (int u = 0; u < 3; ++u) {
		const int v = (u + 1) % 3;
		const int w = (u + 2) % 3;
		const double a = kap[v] * kap[w] / kap[u];
		const double b =
			(sig[v] * kap[w] + sig[w] * kap[v] - a * sig[u]) / kap[u];
		const double c = sig[v] * sig[w] - b * sig[u];
		const double d = sig[u] / kap[u];
		switch (tensor_kind) {
		case 0: // A
			prin[u] = a;
			break;
		case 1: // B
			prin[u] = b;
			break;
		case 2: // C
			prin[u] = c;
			break;
		case 3: // D
			prin[u] = d;
			break;
		case 4: // InvKappa
			prin[u] = 1.0 / kap[u];
			break;
		default:
			prin[u] = 0.0;
			break;
		}
	}

	// Columns of R are principal unit vectors in xyz: ê_r, ê_θ, ê_z (2D) or
	// ê_r, ê_θ, ê_φ (3D). MVP focuses on 2D; 3D uses a simple spherical frame.
	double R[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
	const double inv_r = 1.0 / r;
	if (dim == 2) {
		R[0][0] = dx[0] * inv_r; // ê_r
		R[1][0] = dx[1] * inv_r;
		R[0][1] = -dx[1] * inv_r; // ê_θ
		R[1][1] = dx[0] * inv_r;
		R[2][2] = 1.0; // ê_z
	} else {
		const double rx = dx[0] * inv_r;
		const double ry = dx[1] * inv_r;
		const double rz = dx[2] * inv_r;
		R[0][0] = rx;
		R[1][0] = ry;
		R[2][0] = rz;
		// ê_θ from z-axis cross ê_r, fallback if near poles.
		double tx = -ry;
		double ty = rx;
		double tz = 0.0;
		double tlen = std::sqrt(tx * tx + ty * ty + tz * tz);
		if (tlen < 1e-14) {
			tx = 0.0;
			ty = -rz;
			tz = ry;
			tlen = std::sqrt(tx * tx + ty * ty + tz * tz);
		}
		tx /= tlen;
		ty /= tlen;
		tz /= tlen;
		R[0][1] = tx;
		R[1][1] = ty;
		R[2][1] = tz;
		// ê_φ = ê_r × ê_θ
		R[0][2] = ry * tz - rz * ty;
		R[1][2] = rz * tx - rx * tz;
		R[2][2] = rx * ty - ry * tx;
	}

	// T_xyz = R diag(prin) R^T
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double sum = 0.0;
			for (int k = 0; k < 3; ++k) {
				sum += R[i][k] * prin[k] * R[j][k];
			}
			T_xyz[i][j] = sum;
		}
	}
	return true;
}

void PMLProfileData::buildAttributeMaps(
	mfem::Mesh& mesh, const std::vector<PMLProperties>& regions)
{
	const int max_attr = mesh.attributes.Max();
	is_pml_attr_.assign(max_attr, false);
	attr_region_index_.assign(max_attr, -1);

	for (int ri = 0; ri < static_cast<int>(regions.size()); ++ri) {
		for (const Attribute tag : regions[ri].geom_tags) {
			if (tag <= 0 || tag > max_attr) {
				throw std::runtime_error("PML material tag is out of mesh attribute range.");
			}
			if (is_pml_attr_[tag - 1]) {
				throw std::runtime_error(
					"Overlapping PML material tag assignment for attribute " +
					std::to_string(tag));
			}
			is_pml_attr_[tag - 1] = true;
			attr_region_index_[tag - 1] = ri;
		}
	}

	region_max_depth_.assign(regions.size(), {0.0, 0.0, 0.0});
	radial_.assign(regions.size(), RadialRegionData{});
	global_interfaces_ = {};
}

void PMLProfileData::buildInterfaceData(
	mfem::Mesh& mesh, const std::vector<PMLProperties>& regions)
{
	const int dim = mesh.Dimension();
	(void)regions;

	for (int f = 0; f < mesh.GetNumFaces(); ++f) {
		auto* ft = mesh.GetFaceElementTransformations(f);
		if (!ft || ft->Elem1No < 0 || ft->Elem2No < 0) {
			continue;
		}

		const int attr1 = mesh.GetAttribute(ft->Elem1No);
		const int attr2 = mesh.GetAttribute(ft->Elem2No);
		const bool pml1 = attr1 > 0 && attr1 <= static_cast<int>(is_pml_attr_.size()) &&
		                  is_pml_attr_[attr1 - 1];
		const bool pml2 = attr2 > 0 && attr2 <= static_cast<int>(is_pml_attr_.size()) &&
		                  is_pml_attr_[attr2 - 1];
		if (pml1 == pml2) {
			continue;
		}

		const int pml_attr = pml1 ? attr1 : attr2;
		const int vac_el = pml1 ? ft->Elem2No : ft->Elem1No;
		const int pml_el = pml1 ? ft->Elem1No : ft->Elem2No;

		mfem::Vector vac_center(dim);
		mfem::Vector pml_center(dim);
		mesh.GetElementCenter(vac_el, vac_center);
		mesh.GetElementCenter(pml_el, pml_center);

		mfem::Vector delta(dim);
		delta = pml_center;
		delta -= vac_center;

		mfem::Vector face_center(dim);
		mfem::IntegrationPoint ip;
		ip.Set3(0.5, 0.5, 0.5);
		ft->SetIntPoint(&ip);
		ft->Transform(ip, face_center);

		const int axis = dominantAxis(delta, dim);
		if (axis >= dim) {
			continue;
		}

		auto& iface = global_interfaces_[axis];
		const int sign = (delta(axis) >= 0.0) ? 1 : -1;
		if (sign > 0) {
			// +side: PML lies at larger x_d than the vacuum–PML face.
			if (!iface.set_pos) {
				iface.set_pos = true;
				iface.coord_pos = face_center(axis);
			} else {
				iface.coord_pos = std::min(iface.coord_pos, face_center(axis));
			}
		} else {
			// −side: PML lies at smaller x_d than the vacuum–PML face.
			if (!iface.set_neg) {
				iface.set_neg = true;
				iface.coord_neg = face_center(axis);
			} else {
				iface.coord_neg = std::max(iface.coord_neg, face_center(axis));
			}
		}
	}
}

void PMLProfileData::buildRadialRegionData(
	mfem::Mesh& mesh, const std::vector<PMLProperties>& regions)
{
	const int dim = mesh.Dimension();

	for (int ri = 0; ri < static_cast<int>(regions.size()); ++ri) {
		// True R-stretch and profile-only radial both need center / shell radii.
		if (regions[ri].stretch_mode != PMLStretchMode::Radial &&
		    !regions[ri].uniaxial_radial) {
			continue;
		}

		RadialRegionData& rd = radial_[ri];
		rd.active = true;

		if (regions[ri].radial_center.has_value()) {
			rd.center = *regions[ri].radial_center;
			rd.center_inferred = false;
		} else {
			// Infer center as mean of vacuum–PML interface face centers for this region.
			double sum[3] = {0.0, 0.0, 0.0};
			int n_iface = 0;
			for (int f = 0; f < mesh.GetNumFaces(); ++f) {
				auto* ft = mesh.GetFaceElementTransformations(f);
				if (!ft || ft->Elem1No < 0 || ft->Elem2No < 0) {
					continue;
				}
				const int attr1 = mesh.GetAttribute(ft->Elem1No);
				const int attr2 = mesh.GetAttribute(ft->Elem2No);
				const bool pml1 = attr1 > 0 && attr1 <= static_cast<int>(is_pml_attr_.size()) &&
				                  is_pml_attr_[attr1 - 1];
				const bool pml2 = attr2 > 0 && attr2 <= static_cast<int>(is_pml_attr_.size()) &&
				                  is_pml_attr_[attr2 - 1];
				if (pml1 == pml2) {
					continue;
				}
				const int pml_attr = pml1 ? attr1 : attr2;
				if (attr_region_index_[pml_attr - 1] != ri) {
					continue;
				}
				mfem::Vector face_center(dim);
				mfem::IntegrationPoint ip;
				ip.Set3(0.5, 0.5, 0.5);
				ft->SetIntPoint(&ip);
				ft->Transform(ip, face_center);
				for (int d = 0; d < dim; ++d) {
					sum[d] += face_center(d);
				}
				++n_iface;
			}
			if (n_iface == 0) {
				throw std::runtime_error(
					"PML stretch_mode \"radial\": cannot infer radial_center — "
					"no vacuum–PML interface faces for region " +
					std::to_string(ri) + ". Provide radial_center explicitly.");
			}
			for (int d = 0; d < 3; ++d) {
				rd.center[static_cast<size_t>(d)] =
					(d < dim) ? sum[d] / static_cast<double>(n_iface) : 0.0;
			}
			rd.center_inferred = true;
		}

		rd.r_inner = std::numeric_limits<double>::infinity();
		rd.r_outer = 0.0;
		int n_pml_qp = 0;

		// Bound the shell from PML volume samples. Interface-face centers can sit
		// inward of the true PML volume and under-estimate r_inner, which makes
		// σ(ρ) already large at the geometric vacuum–PML face (looks like a jump
		// even when grading_order ≥ 1) and reflects radial-polarized fields.
		for (int el = 0; el < mesh.GetNE(); ++el) {
			const int attr = mesh.GetAttribute(el);
			if (attr <= 0 || attr > static_cast<int>(is_pml_attr_.size()) ||
			    !is_pml_attr_[attr - 1] || attr_region_index_[attr - 1] != ri) {
				continue;
			}
			mfem::ElementTransformation* T = mesh.GetElementTransformation(el);
			const mfem::IntegrationRule& ir =
				mfem::IntRules.Get(T->GetGeometryType(), fe_order_ + 1);
			for (int iq = 0; iq < ir.GetNPoints(); ++iq) {
				T->SetIntPoint(&ir[iq]);
				mfem::Vector x(dim);
				T->Transform(ir[iq], x);
				const double r = distanceFromCenter(x, rd.center, dim);
				rd.r_inner = std::min(rd.r_inner, r);
				rd.r_outer = std::max(rd.r_outer, r);
				++n_pml_qp;
			}
		}

		if (n_pml_qp == 0 || !std::isfinite(rd.r_inner)) {
			throw std::runtime_error(
				"PML stretch_mode \"radial\": no PML volume quadrature points for region " +
				std::to_string(ri) + ".");
		}
		if (rd.r_outer <= rd.r_inner) {
			throw std::runtime_error(
				"PML stretch_mode \"radial\": non-positive shell thickness for region " +
				std::to_string(ri) + " (r_inner=" + std::to_string(rd.r_inner) +
				", r_outer=" + std::to_string(rd.r_outer) + ").");
		}
	}
}

double PMLProfileData::depthAlongAxis(const mfem::Vector& x, Direction d) const
{
	const auto& iface = global_interfaces_[d];
	double depth = 0.0;
	if (iface.set_pos) {
		depth = std::max(depth, x(d) - iface.coord_pos);
	}
	if (iface.set_neg) {
		depth = std::max(depth, iface.coord_neg - x(d));
	}
	return depth;
}

double PMLProfileData::depthRadial(const mfem::Vector& x, int region_index) const
{
	const RadialRegionData& rd = radial_[region_index];
	const double r = distanceFromCenter(x, rd.center, mesh_dim_);
	return std::max(0.0, r - rd.r_inner);
}

double PMLProfileData::thicknessFor(int region_index, Direction stretch_dir) const
{
	if (regions_[region_index].stretch_mode == PMLStretchMode::Radial &&
	    radial_[region_index].active) {
		return radial_[region_index].thickness();
	}
	return region_max_depth_[region_index][stretch_dir];
}

void PMLProfileData::buildElementProfiles(
	mfem::Mesh& mesh, const std::vector<PMLProperties>& regions, int fe_order)
{
	const int dim = mesh.Dimension();
	const int ne = mesh.GetNE();

	el_to_profile_index_.assign(ne, -1);
	element_profiles_.clear();
	element_profiles_.reserve(ne / 4 + 1);

	for (int el = 0; el < ne; ++el) {
		const int attr = mesh.GetAttribute(el);
		if (attr <= 0 || attr > static_cast<int>(is_pml_attr_.size()) ||
		    !is_pml_attr_[attr - 1]) {
			continue;
		}

		const int region = attr_region_index_[attr - 1];
		const PMLProperties& props = regions[region];

		mfem::ElementTransformation* T = mesh.GetElementTransformation(el);
		const mfem::IntegrationRule& ir =
			mfem::IntRules.Get(T->GetGeometryType(), fe_order + 1);

		PMLElementProfiles ep;
		ep.attribute = attr;
		ep.region_index = region;
		ep.qp_profiles.resize(ir.GetNPoints());

		for (int iq = 0; iq < ir.GetNPoints(); ++iq) {
			T->SetIntPoint(&ir[iq]);
			mfem::Vector x(dim);
			T->Transform(ir[iq], x);

			for (Direction d = X; d <= Z; ++d) {
				if (d >= dim || props.active_axes.count(d) == 0) {
					ep.qp_profiles[iq][d] = PMLDirectionProfiles{};
					continue;
				}

				double rho = 0.0;
				if (props.stretch_mode == PMLStretchMode::Radial && radial_[region].active) {
					rho = depthRadial(x, region);
				} else {
					rho = depthAlongAxis(x, d);
					region_max_depth_[region][d] =
						std::max(region_max_depth_[region][d], rho);
				}
				ep.qp_profiles[iq][d].depth = rho;
			}
			// Uniaxial R: no Cartesian active_axes — store radial depth on slot 0 for diagnostics.
			if (props.uniaxial_radial && radial_[region].active) {
				const double rho = depthRadial(x, region);
				ep.qp_profiles[iq][0].depth = rho;
			}
		}

		element_profiles_.push_back(std::move(ep));
		el_to_profile_index_[el] = static_cast<int>(element_profiles_.size()) - 1;
	}

	for (auto& ep : element_profiles_) {
		const PMLProperties& props = regions[ep.region_index];
		for (int iq = 0; iq < static_cast<int>(ep.qp_profiles.size()); ++iq) {
			if (props.uniaxial_radial) {
				const double L = thicknessFor(ep.region_index, X);
				evaluateStretchProfiles(ep.qp_profiles[iq][0].depth, L, props,
				                        ep.qp_profiles[iq][0]);
				continue;
			}
			for (Direction d = X; d <= Z; ++d) {
				if (d >= dim || props.active_axes.count(d) == 0) {
					continue;
				}
				const double L = thicknessFor(ep.region_index, d);
				evaluateStretchProfiles(ep.qp_profiles[iq][d].depth, L, props,
				                        ep.qp_profiles[iq][d]);
			}
		}
	}
}

void PMLProfileData::printDiagnostics(int rank) const
{
	if (rank != 0 || element_profiles_.empty()) {
		return;
	}

	std::cout << "\n========================================================" << std::endl;
	std::cout << "  VOLUMETRIC PML PROFILE INIT" << std::endl;
	std::cout << "========================================================" << std::endl;
	std::cout << "  PML elements: " << element_profiles_.size() << std::endl;

	for (size_t ri = 0; ri < regions_.size(); ++ri) {
		const auto& props = regions_[ri];
		const char* mode =
			(props.uniaxial_radial)
				? "radial-cylindrical(R)"
				: ((props.stretch_mode == PMLStretchMode::Radial) ? "radial" : "box");
		std::cout << "  Region " << ri << " stretch_mode=" << mode;
		if ((props.stretch_mode == PMLStretchMode::Radial || props.uniaxial_radial) &&
		    ri < radial_.size() && radial_[ri].active) {
			const auto& rd = radial_[ri];
			std::cout << " center=(" << rd.center[0] << ", " << rd.center[1] << ", "
			          << rd.center[2] << ")"
			          << (rd.center_inferred ? " [inferred]" : " [user]")
			          << std::setprecision(8)
			          << " r_inner=" << rd.r_inner << " r_outer=" << rd.r_outer
			          << " L=" << rd.thickness()
			          << std::defaultfloat;
		} else {
			std::cout << " max depth:";
			for (Direction d = X; d <= Z; ++d) {
				std::cout << " axis" << d << "=" << region_max_depth_[ri][d];
			}
			std::cout << " | interfaces:";
			for (Direction d = X; d <= Z; ++d) {
				const auto& iface = global_interfaces_[d];
				if (!iface.set_pos && !iface.set_neg) {
					continue;
				}
				std::cout << " d" << d << "{";
				if (iface.set_neg) {
					std::cout << "−@" << iface.coord_neg;
				}
				if (iface.set_pos && iface.set_neg) {
					std::cout << ",";
				}
				if (iface.set_pos) {
					std::cout << "+@" << iface.coord_pos;
				}
				std::cout << "}";
			}
		}
		std::cout << std::endl;
	}

	double max_iface_sigma = 0.0;
	double max_sigma = 0.0;
	double max_sigma_theta = 0.0;
	double max_iface_sigma_theta = 0.0;
	double max_alpha = 0.0;
	double max_kappa = 1.0;
	bool has_cylindrical_r = false;

	for (const auto& ep : element_profiles_) {
		const PMLProperties& props = regions_[ep.region_index];
		if (props.uniaxial_radial) {
			has_cylindrical_r = true;
			const double L = thicknessFor(ep.region_index, X);
			const double sigma_max = designSigmaMax(L, props);
			const double r_inner =
				(ep.region_index < static_cast<int>(radial_.size()) &&
				 radial_[ep.region_index].active)
					? radial_[ep.region_index].r_inner
					: 0.0;
			for (const auto& qp : ep.qp_profiles) {
				const double rho = qp[0].depth;
				const double r = r_inner + rho;
				const double sig_r = qp[0].sigma;
				const double sig_th = cylindricalSigmaTheta(
					r, rho, L, sigma_max, props.grading_order);
				if (L > 0.0 && rho / L < 0.05) {
					max_iface_sigma = std::max(max_iface_sigma, sig_r);
					max_iface_sigma_theta =
						std::max(max_iface_sigma_theta, sig_th);
				}
				max_sigma = std::max(max_sigma, sig_r);
				max_sigma_theta = std::max(max_sigma_theta, sig_th);
				max_alpha = std::max(max_alpha, qp[0].alpha);
				max_kappa = std::max(max_kappa, qp[0].kappa);
			}
			continue;
		}
		for (const auto& qp : ep.qp_profiles) {
			for (Direction d = X; d <= Z; ++d) {
				const double L = thicknessFor(ep.region_index, d);
				if (L <= 0.0) {
					continue;
				}
				const double rel_depth = qp[d].depth / L;
				if (rel_depth < 0.05) {
					max_iface_sigma = std::max(max_iface_sigma, qp[d].sigma);
				}
				max_sigma = std::max(max_sigma, qp[d].sigma);
				max_alpha = std::max(max_alpha, qp[d].alpha);
				max_kappa = std::max(max_kappa, qp[d].kappa);
			}
		}
	}

	std::cout << std::scientific << std::setprecision(3);
	std::cout << "  Interface-adjacent (depth/L < 0.05): max sigma_r="
	          << max_iface_sigma;
	if (has_cylindrical_r) {
		std::cout << " sigma_theta=" << max_iface_sigma_theta;
	}
	std::cout << std::endl;
	std::cout << "  Global max: sigma_r=" << max_sigma;
	if (has_cylindrical_r) {
		std::cout << " sigma_theta=" << max_sigma_theta;
	}
	std::cout << " kappa=" << max_kappa << " alpha=" << max_alpha << std::endl;
	std::cout << std::defaultfloat;
	std::cout << "========================================================\n" << std::endl;
}

} // namespace maxwell
