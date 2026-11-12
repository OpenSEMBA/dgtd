#pragma once

#include "components/Model.h"
#include "components/ThinWireGeometry.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace maxwell {

// Edelvik in-cell inductance, units of μ. r0 = 1.7 * transverse.
double edelvikNormalizedInductance(double radius, double transverse);
// Cosine tube weight. Zero outside a <= rho <= tube_radius. Cross-section integral is 1.
double edelvikTubeWeight(double rho, double radius, double tube_radius);

std::vector<WirePolyline> buildWirePolylinesFromGmsh(
	const std::string& mesh_path,
	const std::vector<WireSpec>& specs);

// Clip every wire against the serial mesh. The result is the same on every rank.
std::vector<WirePieceGeom> buildWirePieceChain(const Model& model);

class ThinWireCoupler {
public:
	ThinWireCoupler(mfem::ParFiniteElementSpace& fes, const Model& model);

	void commit(double time, double dt);
	void setFeed(const double position[3], double magnitude, double spread, double t0, bool derivative);
	void contribute(const mfem::Vector& in, mfem::Vector& out, int ndofs, double stage_time) const;
	void finalize(const mfem::Vector& eh, int ndofs);

	int lineDofs() const { return static_cast<int>(current_.Size()); }
	bool runsOnDevice() const { return device_path_; }
	int testingOwnedCount() const;
	double testingPieceLength(int piece) const;
	double testingPieceInductance(int piece) const;
	double testingCurrentValue(int dof) const;
	double testingCurrentNorm() const;
	double testingStageCurrentValue(int dof) const;
	double testingStageChargeValue(int dof) const;
	void testingSetUniformCurrent(double value);
	void testingSetOpenHalfWave();
	void testingSetElectricReturnScale(double scale);
	int testingFeedCount() const { return feed_count_; }
	void testingContributeOnHost(const mfem::Vector& in, mfem::Vector& out, int ndofs, double stage_time) const;
	bool sample(const std::vector<double>& point, double& current, double& charge) const;
	void writeParaView(const std::string& directory, int cycle, double time) const;

#ifdef SEMBA_DGTD_ENABLE_CUDA
	// Public so nvcc can emit the extended __device__ lambdas inside it.
	void advanceOnDevice(double tau, const mfem::Vector& eh, mfem::Vector* out, int ndofs,
	                     mfem::Vector& line_current, mfem::Vector& line_charge,
	                     bool copy_from_committed, double sample_time) const;
#endif

private:
	struct TubeBlock {
		mfem::Array<int> vol_dofs;
		mfem::DenseMatrix es_x, es_y, es_z;
		mfem::DenseMatrix j_x, j_y, j_z;
	};

	struct Piece {
		int host = -1;
		bool owned = false;
		double length = 0.0;
		double tangent[3] = {0, 0, 0};
		double inductance = 1.0;
		double extra_inductance = 0.0;
		double resistance = 0.0;
		double v2 = 1.0;
		double radius = 0.0;
		double tube_radius = 0.0;
		int left_neighbor = -1;
		int right_neighbor = -1;
		bool left_pec = false;
		bool right_pec = false;
		int wire = 0;
		double a[3] = {0, 0, 0};
		double b[3] = {0, 0, 0};
		int offset = 0;
		mfem::DenseMatrix mass_inv;
		mfem::DenseMatrix stiffness;
		std::vector<TubeBlock> blocks;
		mfem::Vector shape_left;
		mfem::Vector shape_right;
	};

	static constexpr int kMaxLineDofs = 16;

	void rhs(const mfem::Vector& current, const mfem::Vector& charge,
	         const mfem::Vector* drive,
	         mfem::Vector& d_current, mfem::Vector& d_charge) const;
	void electricDrive(const mfem::Vector* eh, int ndofs, mfem::Vector& drive) const;
	void advance(double tau, const mfem::Vector* eh, int ndofs,
	             mfem::Vector& current, mfem::Vector& charge, double sample_time) const;
	double feedElectric(double time) const;
	void addCurrentToElectric(const mfem::Vector& current, mfem::Vector& out, int ndofs) const;
	void contributeHost(const mfem::Vector& in, mfem::Vector& out, int ndofs, double stage_time) const;
	void packDevice();
	void buildExportMeshes();
	void syncExportFields() const;

	struct WireExport {
		std::vector<int> piece_ids;
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<mfem::DG_FECollection> collection;
		std::unique_ptr<mfem::FiniteElementSpace> space;
		std::unique_ptr<mfem::GridFunction> current;
		std::unique_ptr<mfem::GridFunction> charge;
		std::unique_ptr<mfem::ParaViewDataCollection> paraview;
	};

	const Model& model_;
	MPI_Comm comm_ = MPI_COMM_WORLD;
	int order_ = 1;
	int nvol_ = 0;
	bool device_path_ = false;
	std::vector<Piece> pieces_;
	int nblocks_ = 0;
	mfem::Array<int> block_piece_;
	mfem::Array<int> vol_dofs_;
	mfem::Array<int> offset_;
	mfem::Array<int> left_neighbor_;
	mfem::Array<int> right_neighbor_;
	mfem::Array<int> left_pec_;
	mfem::Array<int> right_pec_;
	mfem::Vector es_;
	mfem::Vector j_mat_;
	mfem::Vector mass_inv_;
	mfem::Vector stiffness_;
	mfem::Vector shape_left_;
	mfem::Vector shape_right_;
	mfem::Vector v2_;
	mfem::Vector inductance_;
	mfem::Vector resistance_;
	mutable mfem::Vector drive_;
	mfem::Vector current_;
	mfem::Vector charge_;
	mutable mfem::Vector stage_current_;
	mutable mfem::Vector stage_charge_;
	mutable mfem::Vector d_current_;
	mutable mfem::Vector d_charge_;
	double step_time_ = 0.0;
	double step_dt_ = 0.0;
	double substep_ = 1.0;
	// Production return is the full transpose. Tests set this to zero to drop it.
	double electric_return_scale_ = 1.0;
	static constexpr int kMaxFeeds = 2;
	int feed_count_ = 0;
	int feed_piece_[kMaxFeeds] = {-1, -1};
	double feed_per_length_[kMaxFeeds] = {0.0, 0.0};
	double feed_magnitude_ = 0.0;
	double feed_spread_ = 1.0;
	double feed_t0_ = 0.0;
	bool feed_derivative_ = false;
	int basis_type_ = mfem::BasisType::GaussLobatto;
	mutable std::vector<WireExport> exports_;
};

}
