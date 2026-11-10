#include <gtest/gtest.h>

#include "wires/ThinWire.h"

#include <cmath>
#include <fstream>
#include <filesystem>
#include <iterator>
#include <mpi.h>
#include <vector>

using namespace maxwell;

TEST(ThinWire, edelvikInductanceAndTubeWeight)
{
	const double coarse = edelvikNormalizedInductance(1e-3, 0.2);
	const double fine = edelvikNormalizedInductance(1e-4, 0.2);
	const double pi = std::acos(-1.0);
	const double tube = 1.7 * 0.2;
	const double expected = std::log((tube + 1e-3) / (2.0 * 1e-3)) / (2.0 * pi);
	EXPECT_NEAR(coarse, expected, 1e-12);
	EXPECT_GT(fine, coarse);

	const double radius = 0.1;
	const double tube_radius = 1.7;
	const int n = 4000;
	const double h = (tube_radius - radius) / static_cast<double>(n);
	double integral = 0.0;
	for (int i = 0; i < n; ++i) {
		const double rho = radius + (static_cast<double>(i) + 0.5) * h;
		integral += edelvikTubeWeight(rho, radius, tube_radius) * 2.0 * pi * rho * h;
	}
	EXPECT_NEAR(integral, 1.0, 1e-4);
	EXPECT_DOUBLE_EQ(0.0, edelvikTubeWeight(0.5 * radius, radius, tube_radius));
	EXPECT_DOUBLE_EQ(0.0, edelvikTubeWeight(tube_radius + 0.01, radius, tube_radius));
}

TEST(ThinWire, uniformTangentialFieldDrivesTheLineAtThatField)
{
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "Tube quadrature check runs on one rank.";
	}
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
		12, 12, 12, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
	Model model(mesh);
	WireSpec spec;
	spec.radius = 1e-3;
	spec.transverse_step = 0.08;
	WirePolyline curve;
	curve.points = {{{0.4, 0.5, 0.5}}, {{0.6, 0.5, 0.5}}};
	model.setWires({spec}, {curve});

	mfem::DG_FECollection collection(2, 3, mfem::BasisType::GaussLobatto);
	mfem::ParFiniteElementSpace space(&model.getMesh(), &collection);
	ThinWireCoupler wire(space, model);

	const int ndofs = space.GetNDofs();
	mfem::Vector field(6 * ndofs);
	field = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		field(i) = 1.0;
	}
	const double dt = 1e-8;
	wire.commit(0.0, dt);
	wire.finalize(field, ndofs);

	const double inductance = wire.testingPieceInductance(0);
	const double expected = dt / inductance;
	for (int i = 0; i < wire.lineDofs(); ++i) {
		EXPECT_NEAR(wire.testingCurrentValue(i), expected, 1e-3 * expected);
	}
}

TEST(ThinWire, readsElementaryTagBeforeVolumeElements)
{
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "The tag fixture writes one mesh file.";
	}
	const std::string path = "thin_wire_elementary_tag.msh";
	{
		std::ofstream out(path);
		out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
		out << "$Nodes\n2\n1 0 0 0\n2 1 0 0\n$EndNodes\n";
		out << "$Elements\n2\n";
		out << "1 1 2 9 4 1 2\n";
		out << "2 4 2 1 1 1 2 3 4\n";
		out << "$EndElements\n";
	}
	WireSpec spec;
	spec.tags = {4};
	spec.radius = 1e-4;
	const auto curves = buildWirePolylinesFromGmsh(path, {spec});
	ASSERT_EQ(1u, curves.size());
	EXPECT_EQ(std::vector<int>({4}), curves[0].tags);
	ASSERT_EQ(2u, curves[0].points.size());
	EXPECT_DOUBLE_EQ(0.0, curves[0].points[0][0]);
	EXPECT_DOUBLE_EQ(1.0, curves[0].points[1][0]);
	std::remove(path.c_str());
}

TEST(ThinWire, openWireCurrentFeedsTheElectricResidual)
{
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "Thin wires are single rank.";
	}
	mfem::Mesh mesh(3, 4, 1, 4);
	mesh.AddVertex(0.0, 0.0, 0.0);
	mesh.AddVertex(1.0, 0.0, 0.0);
	mesh.AddVertex(0.0, 1.0, 0.0);
	mesh.AddVertex(0.0, 0.0, 1.0);
	mesh.AddElement(new mfem::Tetrahedron(0, 1, 2, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 2, 1));
	mesh.AddBdrElement(new mfem::Triangle(0, 1, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 3, 2));
	mesh.AddBdrElement(new mfem::Triangle(1, 2, 3));
	mesh.FinalizeTopology();

	Model model(mesh);
	WireSpec spec;
	spec.radius = 1e-4;
	WirePolyline curve;
	curve.points = {{{0.05, 0.05, 0.05}}, {{0.15, 0.05, 0.05}}};
	model.setWires({spec}, {curve});

	mfem::DG_FECollection collection(1, 3, mfem::BasisType::GaussLobatto);
	mfem::ParFiniteElementSpace space(&model.getMesh(), &collection);
	ThinWireCoupler wire(space, model);
	wire.testingSetUniformCurrent(1.0);
	double sampled_current = 0.0;
	double sampled_charge = 0.0;
	EXPECT_TRUE(wire.sample({0.1, 0.05, 0.05}, sampled_current, sampled_charge));
	EXPECT_NEAR(sampled_current, 1.0, 1e-8);
	EXPECT_NEAR(sampled_charge, 0.0, 1e-8);
	EXPECT_FALSE(wire.sample({0.1, 0.2, 0.05}, sampled_current, sampled_charge));

	const int ndofs = space.GetNDofs();
	mfem::Vector in(6 * ndofs);
	mfem::Vector out(6 * ndofs);
	in = 0.0;
	out = 0.0;
	wire.commit(0.0, 0.0);
	wire.contribute(in, out, ndofs, 0.0);
	out.HostRead();
	const std::string directory = "thin_wire_paraview_test.msh";
	wire.writeParaView(directory, 0, 0.0);
	const std::string pvd =
		"thin_wire_paraview_test.msh/thin_wire_paraview_test_wire00.msh/thin_wire_paraview_test_wire00.msh.pvd";
	std::ifstream pvd_file(pvd);
	std::string contents((std::istreambuf_iterator<char>(pvd_file)), std::istreambuf_iterator<char>());
	EXPECT_NE(std::string::npos, contents.find("Cycle"));
	std::filesystem::remove_all("thin_wire_paraview_test.msh");

	double electric = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		electric += out(i);
	}
	EXPECT_LT(electric, 0.0);

	wire.commit(0.0, 1e-4);
	wire.finalize(in, ndofs);
	EXPECT_LT(wire.testingCurrentNorm(), std::sqrt(static_cast<double>(wire.lineDofs())));
}

TEST(ThinWire, partitionedWireKeepsTheEndpointFlux)
{
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (mfem::Device::Allows(mfem::Backend::CUDA)) {
		GTEST_SKIP() << "CUDA thin wires are single rank.";
	}
#endif
	if (mfem::Mpi::WorldSize() != 2) {
		GTEST_SKIP() << "Requires two MPI ranks.";
	}
	mfem::Mesh mesh(3, 5, 2, 0);
	mesh.AddVertex(0.0, 0.0, 0.0);
	mesh.AddVertex(1.0, 0.0, 0.0);
	mesh.AddVertex(0.0, 1.0, 0.0);
	mesh.AddVertex(0.0, 0.0, 1.0);
	mesh.AddVertex(1.0, 1.0, 1.0);
	mesh.AddElement(new mfem::Tetrahedron(0, 1, 2, 3));
	mesh.AddElement(new mfem::Tetrahedron(1, 2, 3, 4));
	mesh.FinalizeTopology();
	mesh.GenerateBoundaryElements();

	int partition[2] = {0, 1};
	Model model(mesh, GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo{}, partition);
	WireSpec spec;
	spec.radius = 1e-4;
	WirePolyline curve;
	curve.points = {{{0.15, 0.15, 0.15}}, {{0.7, 0.7, 0.7}}};
	model.setWires({spec}, {curve});

	mfem::DG_FECollection collection(1, 3, mfem::BasisType::GaussLobatto);
	mfem::ParFiniteElementSpace space(&model.getMesh(), &collection);
	ThinWireCoupler wire(space, model);
	EXPECT_EQ(4, wire.lineDofs());

	int owned = wire.testingOwnedCount();
	int owned_sum = 0;
	MPI_Allreduce(&owned, &owned_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	EXPECT_EQ(1, owned);
	EXPECT_EQ(2, owned_sum);

	wire.testingSetUniformCurrent(1.0);
	const int ndofs = space.GetNDofs();
	mfem::Vector field(6 * ndofs);
	field = 0.0;
	wire.commit(0.0, 1e-3);
	wire.finalize(field, ndofs);

	const int n = wire.lineDofs();
	std::vector<double> local(static_cast<std::size_t>(n));
	std::vector<double> summed(static_cast<std::size_t>(n));
	for (int i = 0; i < n; ++i) {
		local[static_cast<std::size_t>(i)] = wire.testingCurrentValue(i);
	}
	MPI_Allreduce(local.data(), summed.data(), n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for (int i = 0; i < n; ++i) {
		EXPECT_NEAR(summed[static_cast<std::size_t>(i)], 2.0 * local[static_cast<std::size_t>(i)], 1e-12);
	}
	EXPECT_GT(std::abs(local[0] - 1.0) + std::abs(local[1] - 1.0), 1e-8);
	EXPECT_GT(std::abs(local[2] - 1.0) + std::abs(local[3] - 1.0), 1e-8);

	double sampled_current = 0.0;
	double sampled_charge = 0.0;
	EXPECT_TRUE(wire.sample({0.7, 0.7, 0.7}, sampled_current, sampled_charge));
}

TEST(ThinWire, openHalfWaveFlipsAfterOneLength)
{
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "Line-speed check runs on one rank.";
	}
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
		8, 2, 2, mfem::Element::TETRAHEDRON, 1.6, 0.5, 0.5);
	Model model(mesh);
	WireSpec spec;
	spec.radius = 1e-4;
	WirePolyline curve;
	curve.points = {{{0.3, 0.2, 0.3}}, {{1.3, 0.2, 0.3}}};
	model.setWires({spec}, {curve});

	mfem::DG_FECollection collection(2, 3, mfem::BasisType::GaussLobatto);
	mfem::ParFiniteElementSpace space(&model.getMesh(), &collection);
	ThinWireCoupler wire(space, model);
	wire.testingSetOpenHalfWave();

	double current = 0.0;
	double charge = 0.0;
	ASSERT_TRUE(wire.sample({0.8, 0.2, 0.3}, current, charge));
	const double initial = current;
	ASSERT_GT(initial, 0.5);

	const int ndofs = space.GetNDofs();
	mfem::Vector field(6 * ndofs);
	field = 0.0;
	const double length = 1.0;
	const double dt = 0.01;
	double most_negative = initial;
	double t_negative = 0.0;
	for (int step = 0; step < 250; ++step) {
		const double t = step * dt;
		wire.commit(t, dt);
		wire.finalize(field, ndofs);
		ASSERT_TRUE(wire.sample({0.8, 0.2, 0.3}, current, charge));
		if (current < most_negative) {
			most_negative = current;
			t_negative = t + dt;
		}
	}
	EXPECT_LT(most_negative, -0.5 * initial);
	EXPECT_GT(t_negative, 0.7 * length);
	EXPECT_LT(t_negative, 1.4 * length);
}

namespace {

struct PulseSample {
	double time = 0.0;
	double current = 0.0;
};

std::vector<PulseSample> driveOpenWire(bool electric_return)
{
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
		8, 2, 2, mfem::Element::TETRAHEDRON, 1.6, 0.5, 0.5);
	Model model(mesh);
	WireSpec spec;
	spec.radius = 1e-3;
	WirePolyline curve;
	curve.points = {{{0.3, 0.2, 0.3}}, {{1.3, 0.2, 0.3}}};
	model.setWires({spec}, {curve});

	mfem::DG_FECollection collection(2, 3, mfem::BasisType::GaussLobatto);
	mfem::ParFiniteElementSpace space(&model.getMesh(), &collection);
	ThinWireCoupler wire(space, model);
	const double feed[3] = {0.8, 0.2, 0.3};
	const double spread = 0.05;
	const double t0 = 0.2;
	wire.setFeed(feed, 1.0, spread, t0, false);
	wire.testingSetElectricReturnScale(electric_return ? 1.0 : 0.0);

	const int ndofs = space.GetNDofs();
	mfem::Vector field(6 * ndofs);
	mfem::Vector rate(6 * ndofs);
	field = 0.0;
	const double dt = 0.002;
	const int steps = 800;
	std::vector<PulseSample> samples;
	samples.reserve(static_cast<std::size_t>(steps + 1));
	double current = 0.0;
	double charge = 0.0;
	wire.sample(std::vector<double>(feed, feed + 3), current, charge);
	samples.push_back({0.0, current});
	for (int step = 0; step < steps; ++step) {
		const double t = step * dt;
		wire.commit(t, dt);
		rate = 0.0;
		wire.testingContributeOnHost(field, rate, ndofs, t + dt);
		wire.finalize(field, ndofs);
		if (electric_return) {
			field.Add(dt, rate);
		}
		wire.sample(std::vector<double>(feed, feed + 3), current, charge);
		samples.push_back({t + dt, current});
	}
	return samples;
}

// A returning peak is a local extremum after the drive pulse, near one round trip.
bool echoNearRoundTrip(const std::vector<PulseSample>& samples, double& echo_time, double& echo_current)
{
	const double t0 = 0.2;
	const double spread = 0.05;
	const double pulse_end = t0 + 4.0 * spread;
	double drive_peak = 0.0;
	for (const PulseSample& sample : samples) {
		if (sample.time <= pulse_end) {
			drive_peak = std::max(drive_peak, std::abs(sample.current));
		}
	}
	echo_time = -1.0;
	echo_current = 0.0;
	if (!(drive_peak > 0.0)) {
		return false;
	}
	for (std::size_t i = 1; i + 1 < samples.size(); ++i) {
		const double t = samples[i].time;
		if (t < t0 + 0.7 || t > t0 + 1.4) {
			continue;
		}
		const double a = std::abs(samples[i - 1].current);
		const double b = std::abs(samples[i].current);
		const double c = std::abs(samples[i + 1].current);
		if (b >= a && b >= c && b > 0.2 * drive_peak && b > std::abs(echo_current)) {
			echo_time = t;
			echo_current = samples[i].current;
		}
	}
	return echo_time > 0.0;
}

}  // namespace

TEST(ThinWire, shortPulseEchoWithAndWithoutElectricReturn)
{
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "Pulse echo runs on one rank.";
	}
	const std::vector<PulseSample> coupled = driveOpenWire(true);
	const std::vector<PulseSample> bare = driveOpenWire(false);
	double coupled_time = 0.0;
	double coupled_current = 0.0;
	double bare_time = 0.0;
	double bare_current = 0.0;
	const bool coupled_echo = echoNearRoundTrip(coupled, coupled_time, coupled_current);
	const bool bare_echo = echoNearRoundTrip(bare, bare_time, bare_current);
	std::cout << "pulse echo coupled=" << coupled_echo
	          << " t=" << coupled_time << " I=" << coupled_current
	          << " bare=" << bare_echo
	          << " t=" << bare_time << " I=" << bare_current << "\n";
	for (std::size_t i = 0; i < bare.size(); i += 25) {
		std::cout << "t=" << bare[i].time
		          << " bare=" << bare[i].current
		          << " coupled=" << coupled[i].current << "\n";
	}
	EXPECT_TRUE(bare_echo);
	EXPECT_TRUE(coupled_echo);
}

TEST(ThinWire, pmlHostAllowedWhenKappaIsOne)
{
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "PML host check runs on one rank.";
	}
	mfem::Mesh mesh(3, 4, 1, 4);
	mesh.AddVertex(0.0, 0.0, 0.0);
	mesh.AddVertex(1.0, 0.0, 0.0);
	mesh.AddVertex(0.0, 1.0, 0.0);
	mesh.AddVertex(0.0, 0.0, 1.0);
	mesh.AddElement(new mfem::Tetrahedron(0, 1, 2, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 2, 1));
	mesh.AddBdrElement(new mfem::Triangle(0, 1, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 3, 2));
	mesh.AddBdrElement(new mfem::Triangle(1, 2, 3));
	mesh.FinalizeTopology();
	Model model(mesh);
	WireSpec spec;
	spec.radius = 1e-4;
	WirePolyline curve;
	curve.points = {{{0.05, 0.05, 0.05}}, {{0.15, 0.05, 0.05}}};
	model.setWires({spec}, {curve});

	PMLProperties pml;
	pml.geom_tags = {1};
	pml.kappa_max = 1.0;
	pml.active_axes = {X};
	model.setPMLProperties({pml});
	EXPECT_NO_THROW(buildWirePieceChain(model));

	pml.kappa_max = 2.0;
	model.setPMLProperties({pml});
	EXPECT_THROW(buildWirePieceChain(model), std::runtime_error);
}

#ifdef SEMBA_DGTD_ENABLE_CUDA
TEST(ThinWire, deviceContributeMatchesHost)
{
	if (!mfem::Device::Allows(mfem::Backend::CUDA)) {
		GTEST_SKIP() << "CUDA backend not active.";
	}
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "CUDA thin wires are single rank.";
	}

	mfem::Mesh mesh(3, 4, 1, 4);
	mesh.AddVertex(0.0, 0.0, 0.0);
	mesh.AddVertex(1.0, 0.0, 0.0);
	mesh.AddVertex(0.0, 1.0, 0.0);
	mesh.AddVertex(0.0, 0.0, 1.0);
	mesh.AddElement(new mfem::Tetrahedron(0, 1, 2, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 2, 1));
	mesh.AddBdrElement(new mfem::Triangle(0, 1, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 3, 2));
	mesh.AddBdrElement(new mfem::Triangle(1, 2, 3));
	mesh.FinalizeTopology();

	Model model(mesh);
	WireSpec spec;
	spec.radius = 1e-4;
	WirePolyline curve;
	curve.points = {{{0.05, 0.05, 0.05}}, {{0.15, 0.05, 0.05}}};
	model.setWires({spec}, {curve});

	mfem::DG_FECollection collection(1, 3, mfem::BasisType::GaussLobatto);
	mfem::ParFiniteElementSpace space(&model.getMesh(), &collection);
	ThinWireCoupler wire(space, model);
	ASSERT_TRUE(wire.runsOnDevice());
	ASSERT_GT(wire.lineDofs(), 0);

	const int ndofs = space.GetNDofs();
	mfem::Vector field(6 * ndofs);
	field.UseDevice(true);
	double* host_field = field.HostWrite();
	for (int i = 0; i < field.Size(); ++i) {
		host_field[i] = (i < ndofs) ? 1.0 : 0.0;
	}

	wire.testingSetUniformCurrent(1.0);
	wire.commit(0.0, 1e-3);

	mfem::Vector out_host(6 * ndofs);
	out_host = 0.0;
	wire.testingContributeOnHost(field, out_host, ndofs, 5e-4);
	const int nline = wire.lineDofs();
	std::vector<double> host_current(static_cast<std::size_t>(nline));
	std::vector<double> host_charge(static_cast<std::size_t>(nline));
	for (int i = 0; i < nline; ++i) {
		host_current[static_cast<std::size_t>(i)] = wire.testingStageCurrentValue(i);
		host_charge[static_cast<std::size_t>(i)] = wire.testingStageChargeValue(i);
	}
	out_host.HostRead();
	std::vector<double> host_residual(static_cast<std::size_t>(out_host.Size()));
	for (int i = 0; i < out_host.Size(); ++i) {
		host_residual[static_cast<std::size_t>(i)] = out_host[i];
	}

	mfem::Vector out_device(6 * ndofs);
	out_device.UseDevice(true);
	out_device = 0.0;
	wire.contribute(field, out_device, ndofs, 5e-4);
	for (int i = 0; i < nline; ++i) {
		EXPECT_NEAR(wire.testingStageCurrentValue(i), host_current[static_cast<std::size_t>(i)], 1e-10);
		EXPECT_NEAR(wire.testingStageChargeValue(i), host_charge[static_cast<std::size_t>(i)], 1e-10);
	}
	const double* device_residual = out_device.HostRead();
	for (int i = 0; i < out_device.Size(); ++i) {
		EXPECT_NEAR(device_residual[i], host_residual[static_cast<std::size_t>(i)], 1e-10);
	}
}
#endif

TEST(ThinWire, voltageStepsTheFeedCurrent)
{
	if (mfem::Mpi::WorldSize() != 1) {
		GTEST_SKIP() << "Thin wires are single rank.";
	}
	mfem::Mesh mesh(3, 4, 1, 4);
	mesh.AddVertex(0.0, 0.0, 0.0);
	mesh.AddVertex(1.0, 0.0, 0.0);
	mesh.AddVertex(0.0, 1.0, 0.0);
	mesh.AddVertex(0.0, 0.0, 1.0);
	mesh.AddElement(new mfem::Tetrahedron(0, 1, 2, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 2, 1));
	mesh.AddBdrElement(new mfem::Triangle(0, 1, 3));
	mesh.AddBdrElement(new mfem::Triangle(0, 3, 2));
	mesh.AddBdrElement(new mfem::Triangle(1, 2, 3));
	mesh.FinalizeTopology();

	Model model(mesh);
	WireSpec spec;
	spec.radius = 1e-4;
	WirePolyline curve;
	curve.points = {{{0.05, 0.05, 0.05}}, {{0.15, 0.05, 0.05}}};
	model.setWires({spec}, {curve});

	mfem::DG_FECollection collection(1, 3, mfem::BasisType::GaussLobatto);
	mfem::ParFiniteElementSpace space(&model.getMesh(), &collection);
	ThinWireCoupler wire(space, model);
	const double off_wire[3] = {0.50, 0.50, 0.50};
	EXPECT_THROW(wire.setFeed(off_wire, 1.0, 1.0, 0.0, false), std::runtime_error);
	const double on_wire[3] = {0.10, 0.05, 0.05};
	wire.setFeed(on_wire, 2.0, 1.0e6, 0.0, false);

	const double dt = 1.0e-6;
	const int ndofs = space.GetNDofs();
	mfem::Vector field(6 * ndofs);
	field = 0.0;
	wire.commit(0.0, dt);
	wire.finalize(field, ndofs);

	const double length = wire.testingPieceLength(0);
	const double inductance = wire.testingPieceInductance(0);
	const double voltage = 2.0 * std::exp(-0.5 * std::pow(dt / 1.0e6, 2));
	const double expected = dt * (voltage / length) / inductance;
	ASSERT_EQ(2, wire.lineDofs());
	EXPECT_NEAR(wire.testingCurrentValue(0), expected, 1.0e-12);
	EXPECT_NEAR(wire.testingCurrentValue(1), expected, 1.0e-12);
	double sampled_current = 0.0;
	double sampled_charge = 0.0;
	EXPECT_TRUE(wire.sample({0.10, 0.05, 0.05}, sampled_current, sampled_charge));
	EXPECT_NEAR(sampled_current, expected, 1.0e-12);
	EXPECT_NEAR(sampled_charge, 0.0, 1.0e-12);
}
