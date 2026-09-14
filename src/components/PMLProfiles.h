#pragma once

#include "PMLProperties.h"
#include "Types.h"

#include <algorithm>
#include <array>
#include <vector>

namespace maxwell {

/// Per-quadrature-point stretch coefficients for one active direction.
struct PMLDirectionProfiles {
	double depth = 0.0;
	double sigma = 0.0;
	double kappa = 1.0;
	double alpha = 0.0;
};

/// Profiles at all QPs of one PML element (active directions only in map).
struct PMLElementProfiles {
	Attribute attribute = 0;
	int region_index = -1;
	std::vector<std::array<PMLDirectionProfiles, 3>> qp_profiles;
};

class PMLProfileData {
public:
	PMLProfileData() = default;
	PMLProfileData(mfem::Mesh& mesh, const std::vector<PMLProperties>& regions,
	               int fe_order = 2);

	bool empty() const { return element_profiles_.empty(); }

	const PMLElementProfiles* getElementProfiles(int el) const;

	/// Profile for element/quadrature point (nullptr if not PML or inactive direction).
	const PMLDirectionProfiles* getDirectionProfileAtIP(
		int el, const mfem::IntegrationPoint& ip, Direction stretch_dir) const;

	/// Evaluate stretch profiles at a quadrature point (zero outside PML).
	/// Uses T.Attribute + physical x (MPI-safe); does not use T.ElementNo.
	/// Box: planar depth along stretch_dir. Radial: ρ=max(0,‖x−c‖−r_in) for all axes.
	void evaluateAtTransform(mfem::ElementTransformation& T,
	                         const mfem::IntegrationPoint& ip,
	                         Direction stretch_dir,
	                         PMLDirectionProfiles& out) const;

	int feOrder() const { return fe_order_; }

	void printDiagnostics(int rank) const;

private:
	std::vector<bool> is_pml_attr_;
	std::vector<int> attr_region_index_;
	std::vector<PMLElementProfiles> element_profiles_;
	std::vector<int> el_to_profile_index_;

	void buildAttributeMaps(mfem::Mesh& mesh, const std::vector<PMLProperties>& regions);
	void buildInterfaceData(mfem::Mesh& mesh, const std::vector<PMLProperties>& regions);
	void buildRadialRegionData(mfem::Mesh& mesh, const std::vector<PMLProperties>& regions);
	void buildElementProfiles(mfem::Mesh& mesh, const std::vector<PMLProperties>& regions,
	                          int fe_order);

	double depthAlongAxis(const mfem::Vector& x, Direction d) const;
	double depthRadial(const mfem::Vector& x, int region_index) const;
	double thicknessFor(int region_index, Direction stretch_dir) const;

	/// Planar vacuum–PML interfaces for one stretch axis. Box meshes usually have
	/// both a +side and a −side slab (e.g. top and bottom Y-PML); both must be kept.
	struct InterfaceAxisData {
		bool set_pos = false; ///< PML at larger coordinate than the interface
		bool set_neg = false; ///< PML at smaller coordinate than the interface
		double coord_pos = 0.0;
		double coord_neg = 0.0;
	};

	struct RadialRegionData {
		bool active = false;
		bool center_inferred = false;
		std::array<double, 3> center{{0.0, 0.0, 0.0}};
		double r_inner = 0.0;
		double r_outer = 0.0;
		double thickness() const { return std::max(0.0, r_outer - r_inner); }
	};

	std::array<InterfaceAxisData, 3> global_interfaces_;
	std::vector<std::array<double, 3>> region_max_depth_;
	std::vector<RadialRegionData> radial_;
	std::vector<PMLProperties> regions_;
	int fe_order_ = 2;
	int mesh_dim_ = 0;
};

} // namespace maxwell
