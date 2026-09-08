#include "PMLProfiles.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace maxwell {

namespace {

void evaluateStretchProfiles(
	double rho, double L, const PMLProperties& props, PMLDirectionProfiles& out)
{
	out.depth = rho;
	if (L <= 0.0) {
		out.sigma = 0.0;
		return;
	}

	const double xi = std::clamp(rho / L, 0.0, 1.0);
	const int m = props.grading_order;

	double sigma_max = 0.0;
	if (props.target_reflection > 0.0 && props.target_reflection < 1.0) {
		// m=0: constant σ; m>=1: power-law grade. Same σ_max design formula.
		sigma_max = -(static_cast<double>(m) + 1.0) * std::log(props.target_reflection) /
		            (2.0 * L);
	}

	if (m == 0) {
		// Constant conductivity in the PML volume (abrupt at vacuum interface).
		out.sigma = sigma_max;
	} else {
		out.sigma = sigma_max * std::pow(xi, static_cast<double>(m));
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

} // namespace

PMLProfileData::PMLProfileData(
	mfem::Mesh& mesh, const std::vector<PMLProperties>& regions, int fe_order)
	: fe_order_(fe_order)
{
	if (regions.empty()) {
		return;
	}
	regions_ = regions;

	buildAttributeMaps(mesh, regions);
	buildInterfaceData(mesh, regions);
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

	mfem::ElementTransformation* T = nullptr;
	(void)T;
	(void)ip;
	const int nqp = static_cast<int>(ep->qp_profiles.size());
	if (nqp == 0) {
		return nullptr;
	}

	// Match quadrature point index from reference coordinates (same rule as buildElementProfiles).
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

	const double rho = depthAlongAxis(x, stretch_dir);
	const double L = region_max_depth_[region_index][stretch_dir];
	evaluateStretchProfiles(rho, L, props, out);
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

		const int region = attr_region_index_[pml_attr - 1];
		(void)region;

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
		if (!iface.set) {
			iface.set = true;
			iface.coord = face_center(axis);
			iface.sign_into_pml = sign;
		} else if (iface.sign_into_pml == sign) {
			iface.coord = (sign > 0)
			                  ? std::min(iface.coord, face_center(axis))
			                  : std::max(iface.coord, face_center(axis));
		}
	}
}

double PMLProfileData::depthAlongAxis(const mfem::Vector& x, Direction d) const
{
	const auto& iface = global_interfaces_[d];
	if (!iface.set) {
		return 0.0;
	}

	const double raw = iface.sign_into_pml * (x(d) - iface.coord);
	return std::max(0.0, raw);
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

				const double rho = depthAlongAxis(x, d);
				region_max_depth_[region][d] =
					std::max(region_max_depth_[region][d], rho);
				ep.qp_profiles[iq][d].depth = rho;
			}
		}

		element_profiles_.push_back(std::move(ep));
		el_to_profile_index_[el] = static_cast<int>(element_profiles_.size()) - 1;
	}

	for (auto& ep : element_profiles_) {
		const PMLProperties& props = regions[ep.region_index];
		for (int iq = 0; iq < static_cast<int>(ep.qp_profiles.size()); ++iq) {
			for (Direction d = X; d <= Z; ++d) {
				if (d >= dim || props.active_axes.count(d) == 0) {
					continue;
				}
				const double L = region_max_depth_[ep.region_index][d];
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

	for (size_t ri = 0; ri < region_max_depth_.size(); ++ri) {
		std::cout << "  Region " << ri << " max depth:";
		for (Direction d = X; d <= Z; ++d) {
			std::cout << " axis" << d << "=" << region_max_depth_[ri][d];
		}
		std::cout << std::endl;
	}

	double max_iface_sigma = 0.0;
	double max_sigma = 0.0;

	for (const auto& ep : element_profiles_) {
		for (const auto& qp : ep.qp_profiles) {
			for (Direction d = X; d <= Z; ++d) {
				const double L = region_max_depth_[ep.region_index][d];
				if (L <= 0.0) {
					continue;
				}
				const double rel_depth = qp[d].depth / L;
				if (rel_depth < 0.05) {
					max_iface_sigma = std::max(max_iface_sigma, qp[d].sigma);
				}
				max_sigma = std::max(max_sigma, qp[d].sigma);
			}
		}
	}

	std::cout << std::scientific << std::setprecision(3);
	std::cout << "  Interface-adjacent (depth/L < 0.05): max sigma="
	          << max_iface_sigma << std::endl;
	std::cout << "  Global max: sigma=" << max_sigma << std::endl;
	std::cout << std::defaultfloat;
	std::cout << "========================================================\n" << std::endl;
}

} // namespace maxwell
