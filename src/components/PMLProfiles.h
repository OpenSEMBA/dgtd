#pragma once

#include "PMLProperties.h"
#include "Types.h"

#include <vector>

namespace maxwell {

/// Per-quadrature-point stretch coefficients for one active direction.
struct PMLDirectionProfiles {
	double depth = 0.0;
	double sigma = 0.0;
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
	void buildElementProfiles(mfem::Mesh& mesh, const std::vector<PMLProperties>& regions,
	                          int fe_order);

	double depthAlongAxis(const mfem::Vector& x, Direction d) const;

	struct InterfaceAxisData {
		bool set = false;
		double coord = 0.0;
		int sign_into_pml = 1;
	};

	std::array<InterfaceAxisData, 3> global_interfaces_;
	std::vector<std::array<double, 3>> region_max_depth_;
	std::vector<PMLProperties> regions_;
	int fe_order_ = 2;
};

} // namespace maxwell
