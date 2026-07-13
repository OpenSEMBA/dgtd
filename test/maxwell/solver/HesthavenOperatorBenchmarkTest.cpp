#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstring>

#include "TestUtils.h"
#include "components/DGOperatorFactory.h"
#include "components/Model.h"
#include "driver/driver.h"
#include "evolution/GlobalEvolution.h"
#include "evolution/HesthavenEvolution.h"
#include "evolution/Fields.h"
#include "solver/SourcesManager.h"

using namespace mfem;
using namespace maxwell;
using namespace maxwell::driver;

namespace {

std::unique_ptr<Model> makePECModel(Mesh& serial_mesh)
{
	GeomTagToBoundary pecBdr{{1, BdrCond::PEC}};
	return std::make_unique<Model>(
		serial_mesh, GeomTagToMaterialInfo{},
		GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));
}

void fillExtendedInputWithFaceNeighbors(ParFiniteElementSpace& fes, const Vector& x, Vector& x_ext)
{
	const int ndofs = fes.GetNDofs();
	const int nbr = fes.num_face_nbr_dofs;
	const int block = ndofs + nbr;
	std::array<ParGridFunction, 3> e_gf, h_gf;
	for (int d = X; d <= Z; ++d) {
		e_gf[d].SetSpace(&fes);
		h_gf[d].SetSpace(&fes);
		e_gf[d].MakeRef(&fes, const_cast<double*>(x.HostRead()) + d * ndofs);
		h_gf[d].MakeRef(&fes, const_cast<double*>(x.HostRead()) + (3 + d) * ndofs);
	}
	for (int d = X; d <= Z; ++d) {
		e_gf[d].ExchangeFaceNbrData();
		h_gf[d].ExchangeFaceNbrData();
	}
	for (int d = 0; d < 6; ++d) {
		const int field = d % 3;
		const bool is_e = d < 3;
		std::memcpy(x_ext.GetData() + d * block, x.GetData() + d * ndofs,
		            static_cast<size_t>(ndofs) * sizeof(double));
		if (nbr > 0) {
			const ParGridFunction& nbr_gf = is_e ? e_gf[field] : h_gf[field];
			std::memcpy(x_ext.GetData() + d * block + ndofs,
			            nbr_gf.FaceNbrData().HostRead(),
			            static_cast<size_t>(nbr) * sizeof(double));
		}
	}
}

} // namespace

namespace {

double hostL2Distance(const mfem::Vector& a, const mfem::Vector& b)
{
	const double* ah = a.HostRead();
	const double* bh = b.HostRead();
	double dist = 0.0;
	for (int i = 0; i < a.Size(); ++i) {
		const double d = ah[i] - bh[i];
		dist += d * d;
	}
	return std::sqrt(dist);
}

double hostL2Norm(const mfem::Vector& v)
{
	const double* h = v.HostRead();
	double n = 0.0;
	for (int i = 0; i < v.Size(); ++i) {
		n += h[i] * h[i];
	}
	return std::sqrt(n);
}

} // namespace

TEST(HesthavenOperatorBenchmark, Global_vs_Hesthaven_Mult_agreement)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Operator agreement reference uses undivided mesh; run with np=1.";
	}
	const int nx = 4;
	const int ny = 4;
	const int order = 2;

	Mesh serial_mesh = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE, true, 1.0, 1.0);
	auto model = makePECModel(serial_mesh);
	DG_FECollection fec(order, 2, BasisType::GaussLobatto);
	ParFiniteElementSpace fes(&model->getMesh(), &fec);

	Sources sources;
	Fields<ParFiniteElementSpace, ParGridFunction> fields(fes);
	SourcesManager srcmngr(sources, fes, fields);
	EvolutionOptions opts;
	opts.order = order;
	opts.alpha = 1.0;

	Probes probes;
	ProblemDescription pd(*model, probes, sources, opts);
	DGOperatorFactory<ParFiniteElementSpace> dgops(pd, fes);
	auto global_op = dgops.buildGlobalOperator();

	HesthavenEvolution hest(fes, *model, srcmngr, opts);

	const int ndofs = fes.GetNDofs();
	const int nbr = fes.num_face_nbr_dofs;
	const int block = ndofs + nbr;
	Vector x(6 * ndofs), y_hest(6 * ndofs), y_global(6 * ndofs), x_ext(6 * block);
	x = 0.0;
	for (int i = 0; i < x.Size(); ++i) {
		x[i] = 0.1 * std::sin(0.3 * i);
	}
	for (int d = 0; d < 6; ++d) {
		std::memcpy(x_ext.GetData() + d * block, x.GetData() + d * ndofs,
		            static_cast<size_t>(ndofs) * sizeof(double));
	}
	fillExtendedInputWithFaceNeighbors(fes, x, x_ext);

	hest.Mult(x, y_hest);
	global_op->Mult(x_ext, y_global);

	y_hest.HostRead();
	y_global.HostRead();

	const double dist = y_hest.DistanceTo(y_global);
	if (Mpi::WorldRank() == 0) {
		std::cout << "[HesthavenOperatorBenchmark] ||y_hest - y_global||_2 = "
		          << dist << "\n";
	}
	EXPECT_LT(dist, 5e-2);
}

TEST(HesthavenOperatorBenchmark, MPI_partition_Mult_is_finite)
{
	if (Mpi::WorldSize() < 2) {
		GTEST_SKIP() << "Requires MPI world size >= 2.";
	}

	const int nx = 6;
	const int ny = 6;
	const int order = 1;

	Mesh serial_mesh = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE, true, 1.0, 1.0);
	auto model = makePECModel(serial_mesh);
	DG_FECollection fec(order, 2, BasisType::GaussLobatto);
	ParFiniteElementSpace fes(&model->getMesh(), &fec);

	Sources sources;
	Fields<ParFiniteElementSpace, ParGridFunction> fields(fes);
	SourcesManager srcmngr(sources, fes, fields);
	EvolutionOptions opts;
	opts.order = order;

	HesthavenEvolution hest(fes, *model, srcmngr, opts);

	const int ndofs = fes.GetNDofs();
	Vector x(6 * ndofs), y(6 * ndofs);
	x = 1.0;
	hest.Mult(x, y);

	y.HostRead();
	double local_norm = y.Norml2();
	double global_norm = 0.0;
	MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM,
	              model->getMesh().GetComm());
	EXPECT_TRUE(std::isfinite(global_norm));
	EXPECT_GT(global_norm, 0.0);
}

TEST(HesthavenOperatorBenchmark, MPI_partition_matches_global_operator)
{
	if (Mpi::WorldSize() < 2) {
		GTEST_SKIP() << "Requires MPI world size >= 2.";
	}

	const int nx = 4;
	const int ny = 4;
	const int order = 2;

	Mesh serial_mesh = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE, true, 1.0, 1.0);
	auto model = makePECModel(serial_mesh);
	DG_FECollection fec(order, 2, BasisType::GaussLobatto);
	ParFiniteElementSpace fes(&model->getMesh(), &fec);

	Sources sources;
	Fields<ParFiniteElementSpace, ParGridFunction> fields(fes);
	SourcesManager srcmngr(sources, fes, fields);
	EvolutionOptions opts;
	opts.order = order;
	opts.alpha = 1.0;

	Probes probes;
	ProblemDescription pd(*model, probes, sources, opts);
	DGOperatorFactory<ParFiniteElementSpace> dgops(pd, fes);
	auto global_op = dgops.buildGlobalOperator();

	HesthavenEvolution hest(fes, *model, srcmngr, opts);

	const int ndofs = fes.GetNDofs();
	const int nbr = fes.num_face_nbr_dofs;
	const int block = ndofs + nbr;
	Vector x(6 * ndofs), y_hest(6 * ndofs), y_global(6 * ndofs), x_ext(6 * block);
	x = 0.0;
	for (int i = 0; i < x.Size(); ++i) {
		x[i] = 0.1 * std::sin(0.17 * i + Mpi::WorldRank());
	}
	for (int d = 0; d < 6; ++d) {
		std::memcpy(x_ext.GetData() + d * block, x.GetData() + d * ndofs,
		            static_cast<size_t>(ndofs) * sizeof(double));
	}
	fillExtendedInputWithFaceNeighbors(fes, x, x_ext);

	hest.Mult(x, y_hest);
	global_op->Mult(x_ext, y_global);

	y_hest.HostRead();
	y_global.HostRead();

	const double dist = y_hest.DistanceTo(y_global);
	const double norm_global = y_global.Norml2();
	double global_dist = 0.0;
	double global_norm = 0.0;
	MPI_Allreduce(&dist, &global_dist, 1, MPI_DOUBLE, MPI_MAX, model->getMesh().GetComm());
	MPI_Allreduce(&norm_global, &global_norm, 1, MPI_DOUBLE, MPI_MAX, model->getMesh().GetComm());
	if (Mpi::WorldRank() == 0) {
		std::cout << "[HesthavenOperatorBenchmark MPI] max rank ||y_hest - y_global||_2 = "
		          << global_dist << " (rel " << global_dist / std::max(global_norm, 1e-30) << ")\n";
	}
	EXPECT_LT(global_dist / std::max(global_norm, 1e-30), 5e-2);
}

TEST(HesthavenOperatorBenchmark, MPI_GlobalEvolution_matches_global_operator)
{
	if (Mpi::WorldSize() < 2) {
		GTEST_SKIP() << "Requires MPI world size >= 2.";
	}

	const int nx = 4;
	const int ny = 4;
	const int order = 2;

	Mesh serial_mesh = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE, true, 1.0, 1.0);
	auto model = makePECModel(serial_mesh);
	DG_FECollection fec(order, 2, BasisType::GaussLobatto);
	ParFiniteElementSpace fes(&model->getMesh(), &fec);

	Sources sources;
	Fields<ParFiniteElementSpace, ParGridFunction> fields(fes);
	SourcesManager srcmngr(sources, fes, fields);
	EvolutionOptions opts;
	opts.order = order;
	opts.alpha = 1.0;

	Probes probes;
	ProblemDescription pd(*model, probes, sources, opts);
	DGOperatorFactory<ParFiniteElementSpace> dgops(pd, fes);
	auto global_op = dgops.buildGlobalOperator();

	GlobalEvolution global_evol(fes, *model, srcmngr, opts, probes, 1.0);

	const int ndofs = fes.GetNDofs();
	const int nbr = fes.num_face_nbr_dofs;
	const int block = ndofs + nbr;
	Vector x(6 * ndofs), y_evol(6 * ndofs), y_global(6 * ndofs), x_ext(6 * block);
	x = 0.0;
	for (int i = 0; i < x.Size(); ++i) {
		x[i] = 0.1 * std::sin(0.17 * i + Mpi::WorldRank());
	}
	for (int d = 0; d < 6; ++d) {
		std::memcpy(x_ext.GetData() + d * block, x.GetData() + d * ndofs,
		            static_cast<size_t>(ndofs) * sizeof(double));
	}
	fillExtendedInputWithFaceNeighbors(fes, x, x_ext);

	global_evol.Mult(x, y_evol);
	global_op->Mult(x_ext, y_global);

	y_evol.HostRead();
	y_global.HostRead();

	const double dist = y_evol.DistanceTo(y_global);
	const double norm_global = y_global.Norml2();
	double global_dist = 0.0;
	double global_norm = 0.0;
	MPI_Allreduce(&dist, &global_dist, 1, MPI_DOUBLE, MPI_MAX, model->getMesh().GetComm());
	MPI_Allreduce(&norm_global, &global_norm, 1, MPI_DOUBLE, MPI_MAX, model->getMesh().GetComm());
	if (Mpi::WorldRank() == 0) {
		std::cout << "[GlobalEvolution MPI] max rank ||y_evol - y_global||_2 = "
		          << global_dist << " (rel " << global_dist / std::max(global_norm, 1e-30) << ")\n";
	}
	EXPECT_LT(global_dist / std::max(global_norm, 1e-30), 5e-2);
}

#ifdef SEMBA_DGTD_ENABLE_CUDA
TEST(HesthavenOperatorBenchmark, GPU_Mult_matches_CPU_Mult)
{
	if (!Device::Allows(Backend::CUDA)) {
		GTEST_SKIP() << "CUDA backend not active.";
	}
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "GPU agreement reference uses undivided mesh; run with np=1.";
	}

	const int nx = 4;
	const int ny = 4;
	const int order = 2;

	Mesh serial_mesh = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE, true, 1.0, 1.0);
	auto model = makePECModel(serial_mesh);
	DG_FECollection fec(order, 2, BasisType::GaussLobatto);
	ParFiniteElementSpace fes(&model->getMesh(), &fec);

	Sources sources;
	Fields<ParFiniteElementSpace, ParGridFunction> fields(fes);
	SourcesManager srcmngr(sources, fes, fields);
	EvolutionOptions opts;
	opts.order = order;
	opts.alpha = 1.0;

	HesthavenEvolution hest(fes, *model, srcmngr, opts);

	const int ndofs = fes.GetNDofs();
	Vector x(6 * ndofs), y_cpu(6 * ndofs), y_gpu(6 * ndofs);
	y_gpu.UseDevice(true);
	x = 0.0;
	for (int i = 0; i < x.Size(); ++i) {
		x[i] = 0.1 * std::sin(0.31 * i);
	}

	hest.benchmarkMultCPU(x, y_cpu);
	hest.Mult(x, y_gpu);

	y_cpu.HostRead();
	y_gpu.HostRead();
	const double dist = hostL2Distance(y_cpu, y_gpu);
	const double norm_cpu = hostL2Norm(y_cpu);
	if (Mpi::WorldRank() == 0) {
		std::cout << "[Hesthaven GPU] ||y_cpu - y_gpu||_2 = " << dist
		          << " (rel " << dist / std::max(norm_cpu, 1e-30)
		          << ", norm_cpu=" << norm_cpu << ", norm_gpu=" << hostL2Norm(y_gpu) << ")\n";
	}
	EXPECT_LT(dist / std::max(norm_cpu, 1e-30), 1e-10);
}

TEST(HesthavenOperatorBenchmark, GPU_Mult_TFSF_matches_CPU_Mult)
{
	if (!Device::Allows(Backend::CUDA)) {
		GTEST_SKIP() << "CUDA backend not active.";
	}
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "GPU agreement reference uses undivided mesh; run with np=1.";
	}

	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_IntBoundary"));
	case_data["solver_options"]["evolution_operator"] = "hesthaven";
	auto solver = buildSolver(case_data, maxwellCase("2D_TFSF_IntBoundary"), true);
	auto* hest = dynamic_cast<HesthavenEvolution*>(solver.getEvolTDO());
	ASSERT_NE(nullptr, hest);

	const int ndofs = solver.getFES().GetNDofs();
	Vector x(6 * ndofs), y_cpu(6 * ndofs), y_gpu(6 * ndofs);
	x.UseDevice(true);
	y_cpu.UseDevice(false);
	y_gpu.UseDevice(true);
	x = 0.0;
	for (int i = 0; i < x.Size(); ++i) {
		x.HostWrite()[i] = 0.1 * std::sin(0.23 * i);
	}

	hest->SetTime(0.25);
	hest->benchmarkMultCPU(x, y_cpu);
	hest->Mult(x, y_gpu);

	y_cpu.HostRead();
	y_gpu.HostRead();
	const double dist = hostL2Distance(y_cpu, y_gpu);
	const double norm_cpu = hostL2Norm(y_cpu);
	if (Mpi::WorldRank() == 0) {
		std::cout << "[Hesthaven GPU TFSF] ||y_cpu - y_gpu||_2 = " << dist
		          << " (rel " << dist / std::max(norm_cpu, 1e-30) << ")\n";
	}
	EXPECT_LT(dist / std::max(norm_cpu, 1e-30), 1e-9);
}
#endif
