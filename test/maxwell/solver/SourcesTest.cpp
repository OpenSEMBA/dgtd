#include <gtest/gtest.h>

#include "SourceFixtures.h"

#include "TestUtils.h"
#include "driver/driver.h"
#include "solver/Solver.h"

using namespace maxwell;
using namespace maxwell::driver;
using namespace mfem;
using namespace fixtures::sources;

using Solver = maxwell::Solver;

class SourcesTest : public ::testing::Test {
protected:
	
};

TEST_F(SourcesTest, planewave)
{
	Planewave pw{Gaussian{ 0.2, mfem::Vector({-1.5}) }, unitVec(Y), unitVec(X), E};

	const double TOL{ 1e-6 };

	// Follows maximum.
	EXPECT_NEAR(1.0, pw.eval(Source::Position({ -1.5 }), 0.0, E, Y), TOL);
	EXPECT_NEAR(1.0, pw.eval(Source::Position({  0.0 }), 1.5, E, Y), TOL);
	
	// H has same magnitude as E.
	for (auto x{0.0}; x < 3.0; x += 0.1) {
		EXPECT_NEAR(
			pw.eval(Source::Position({ x }), 0.0, E, Y), 
			pw.eval(Source::Position({ x }), 0.0, H, Z), 
			TOL
		);
	}	

}

TEST_F(SourcesTest, DirectTFSFPlanewaveMatchesProjectedReferenceOnCuda)
{
	if (!Device::Allows(Backend::CUDA)) {
		GTEST_SKIP() << "CUDA backend not active.";
	}
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Reference comparison expects an undivided mesh.";
	}

	auto case_data = parseJSONfile(maxwellCase("2D_TFSF"));
	case_data["solver_options"]["evolution_operator"] = "global";
	auto solver = buildSolver(case_data, maxwellCase("2D_TFSF"), true);

	auto& srcmngr = solver.getSourcesManager();
	ASSERT_TRUE(srcmngr.hasDirectEval());

	const double time = 0.25;
	srcmngr.evalTimeVarFieldDirect(time);
	const auto& direct = srcmngr.getCachedTFSFFields();

	Fields<ParFiniteElementSpace, ParGridFunction> ref_fields(solver.getFES());
	SourcesManager ref_srcmngr(srcmngr.sources, solver.getFES(), ref_fields);
	const auto& tfsf_marker =
		solver.getModel().getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn);
	ref_srcmngr.initTFSFPreReqs(solver.getModel().getConstMesh(), tfsf_marker);

	auto projected = ref_srcmngr.evalTimeVarField(time, ref_srcmngr.getGlobalTFSFSpace());
	ref_srcmngr.markDoFSforTForSF(projected, true);
	if (ref_srcmngr.getTFSFSubMesher().getSFSubMesh() != NULL) {
		auto projected_sf = ref_srcmngr.evalTimeVarField(time, ref_srcmngr.getGlobalTFSFSpace());
		ref_srcmngr.markDoFSforTForSF(projected_sf, false);
		for (int f : {E, H}) {
			for (int d : {X, Y, Z}) {
				projected[f][d] -= projected_sf[f][d];
				projected[f][d] *= 0.5;
			}
		}
	}

	const int ndofs = srcmngr.getGlobalTFSFSpace()->GetNDofs();
	const int sample_count = std::min(ndofs, 16);
	const double tol = 1e-11;
	for (int f : {E, H}) {
		for (int d : {X, Y, Z}) {
			const double* direct_data = direct[f][d].HostRead();
			const double* projected_data = projected[f][d].HostRead();
			for (int i = 0; i < sample_count; ++i) {
				EXPECT_NEAR(direct_data[i], projected_data[i], tol)
					<< "Mismatch at dof " << i << " for field " << f << " component " << d;
			}
		}
	}
}
