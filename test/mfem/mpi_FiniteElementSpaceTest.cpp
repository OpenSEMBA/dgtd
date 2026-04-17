#include <gtest/gtest.h>
#include <mfem.hpp>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <vector>

#include "math/PhysicalConstants.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "components/Model.h"

using namespace mfem;

namespace maxwell{

using Rank = int;
using GlobalDoF = int;
using LocalDoF = int;

class MPIFiniteElementSpaceTest : public ::testing::Test {
protected:

	using NodeCoordinate = std::vector<double>;
	using CollocatedNodes = std::map<int, int>;

	std::map<NodeCoordinate, int> buildNodeCoordinateToDofs(
		int elem, const FiniteElementSpace& fes, const GridFunction& nodes)
	{
		std::map<NodeCoordinate, int> localNodesToDof;
		Array<int> localDofs;
		fes.GetElementDofs(elem, localDofs);
		for (int i{ 0 }; i < localDofs.Size(); ++i) {
			int localDof{ localDofs[i] };
			const auto dim{ fes.GetVDim() };
			NodeCoordinate node(dim);
			for (int d{ 0 }; d < dim; ++d) {
				node[d] = nodes[localDof + d];
			}
			localNodesToDof[node] = localDof;
		}
		return localNodesToDof;
	}

	CollocatedNodes getCollocatedNodes(FiniteElementSpace& fes)
	{
		GridFunction nodes{ &fes };
		auto& m{ *fes.GetMesh() };
		m.GetNodes(nodes);

		CollocatedNodes res;
		for (int f{ 0 }; f < m.GetNumFaces(); ++f) {
			const auto faceInfo{ m.GetFaceInformation(f) };
			if (!faceInfo.IsInterior()) {
				continue;
			}
			const auto localNodes{
				buildNodeCoordinateToDofs(faceInfo.element[0].index, fes, nodes) };
			const auto neighNodes{
				buildNodeCoordinateToDofs(faceInfo.element[1].index, fes, nodes) };

			for (const auto& [node, localDof] : localNodes) {
				const auto it{ neighNodes.find(node) };
				if (it != neighNodes.end()) {
					res.emplace(localDof, it->second);
				}
			}
		}
		return res;
	}
};

TEST_F(MPIFiniteElementSpaceTest, GlobalToLocalPartitionElementMaps_np2){

	if (mfem::Mpi::WorldSize() != 2){
		std::cout << "This test must be run with exactly 2 MPI ranks.\n"
				<< "Run it using: mpirun -np 2 <executable>\n";   // or abort if needed
		ASSERT_TRUE(false);
	}
	
	auto mesh = Mesh::MakeCartesian1D(4);
	auto fec = DG_FECollection(1, 1, BasisType::GaussLobatto);
	auto fes = FiniteElementSpace(&mesh, &fec);
	int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
	auto pmesh = ParMesh(MPI_COMM_WORLD, mesh, partitioning);
	auto pfes = ParFiniteElementSpace(&pmesh,&fec);

	auto serial_element_to_center = buildSerialElem2CenterMap(*fes.GetMesh());
	auto part_element_to_center = buildPartitionElem2CenterMap(*pfes.GetParMesh());
	auto g2l_element_map = buildGlobalToPartitionLocalElementMap(serial_element_to_center, part_element_to_center);

	if (Mpi::WorldRank() == 0){
		ASSERT_EQ( 0, g2l_element_map[2]);
		ASSERT_EQ( 1, g2l_element_map[3]);
		std::cout << "Rank 0 OK." << std::endl;
	}
	else if (Mpi::WorldRank() == 1){
		ASSERT_EQ( 0, g2l_element_map[0]);
		ASSERT_EQ( 1, g2l_element_map[1]);
		std::cout << "Rank 1 OK." << std::endl;
	}

}

TEST_F(MPIFiniteElementSpaceTest, GlobalToLocalPartitionDoFMaps_np2)
{
	if (mfem::Mpi::WorldSize() != 2){
		std::cout << "This test must be run with exactly 2 MPI ranks.\n"
				<< "Run it using: mpirun -np 2 <executable>\n";   // or abort if needed
		ASSERT_TRUE(false);
	}
	
	auto mesh = Mesh::MakeCartesian1D(4);
	auto fec = DG_FECollection(1, 1, BasisType::GaussLobatto);
	auto fes = FiniteElementSpace(&mesh, &fec);
	int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
	auto pmesh = ParMesh(MPI_COMM_WORLD, mesh, partitioning);
	auto pfes = ParFiniteElementSpace(&pmesh,&fec);

	pfes.ExchangeFaceNbrData();

	auto serial_element_to_center = buildSerialElem2CenterMap(*fes.GetMesh());
	auto part_element_to_center = buildPartitionElem2CenterMap(*pfes.GetParMesh());
	auto g2l_element_map = buildGlobalToPartitionLocalElementMap(serial_element_to_center, part_element_to_center);

	auto global_dof_table = fes.GetElementToDofTable();
	std::map<GlobalElementId, Array<int>> temp_global_map;
	for (auto e = 0; e < global_dof_table.Size(); e++){
		Array<int> dofs;
		global_dof_table.GetRow(e, dofs);
		temp_global_map[e] = dofs;
	}
	auto local_dof_table = pfes.GetElementToDofTable();
	std::map<LocalElementId, Array<int>> temp_local_map;
	for (auto e = 0; e < local_dof_table.Size(); e++){
		Array<int> dofs;
		local_dof_table.GetRow(e, dofs);
		temp_local_map[e] = dofs;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	std::map<GlobalDoF, LocalDoF> g2l_dof_map;
	for (const auto& [global_el, local_el] : g2l_element_map){
		const auto& global_el_dofs = temp_global_map[global_el];
		const auto& local_el_dofs = temp_local_map[local_el];
		assert(global_el_dofs.Size() == local_el_dofs.Size());
		for (auto d = 0; d < global_el_dofs.Size(); d++){
			g2l_dof_map[global_el_dofs[d]] = local_el_dofs[d];
		}
	}

	if (Mpi::WorldRank() == 0){
		ASSERT_EQ( 0, g2l_dof_map[4]);
		ASSERT_EQ( 1, g2l_dof_map[5]);
		ASSERT_EQ( 2, g2l_dof_map[6]);
		ASSERT_EQ( 3, g2l_dof_map[7]);
		std::cout << "Rank 0 OK." << std::endl;
	}
	else if (Mpi::WorldRank() == 1){
		ASSERT_EQ( 0, g2l_dof_map[0]);
		ASSERT_EQ( 1, g2l_dof_map[1]);
		ASSERT_EQ( 2, g2l_dof_map[2]);
		ASSERT_EQ( 3, g2l_dof_map[3]);
		
		std::cout << "Rank 1 OK." << std::endl;
	}

}

TEST_F(MPIFiniteElementSpaceTest, Scaling2D)
{
	auto mesh = Mesh::MakeCartesian2D(5, 5, Element::Type::QUADRILATERAL, true, 5.0, 5.0);
	auto fec = DG_FECollection(1, 2, BasisType::GaussLobatto);
	auto fes = FiniteElementSpace(&mesh, &fec);
	auto world_size = Mpi::WorldSize();
	int* partitioning = mesh.GeneratePartitioning(world_size);
	auto pmesh = ParMesh(MPI_COMM_WORLD, mesh, partitioning);
	auto pfes = ParFiniteElementSpace(&pmesh,&fec);

	auto serial_element_to_center = buildSerialElem2CenterMap(*fes.GetMesh());
	auto part_element_to_center = buildPartitionElem2CenterMap(*pfes.GetParMesh());
	auto g2l_element_map = buildGlobalToPartitionLocalElementMap(serial_element_to_center, part_element_to_center);

	auto ndof_per_geom = 4; //Order 1 quadrilaterals = 4 dofs;
	std::vector<int> global_dof_indices;
	global_dof_indices.reserve(4 * g2l_element_map.size());
	for (const auto& [global, vals] : g2l_element_map){
		for (auto d = 0; d < ndof_per_geom; d++){
			global_dof_indices.push_back(global * ndof_per_geom + d);
		}
	}

	ConstantCoefficient one(1.0);

	// Serial part

	auto bf = BilinearForm(&fes);
	bf.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
	bf.AddInteriorFaceIntegrator(new maxwell::mfemExtension::MaxwellDGOneNormalJumpIntegrator({maxwell::X}, 1.0));
	bf.Assemble();
	bf.Finalize();

	GridFunction serial_gf(&fes);
	for (auto v = 0; v < serial_gf.Size(); v++){
		serial_gf[v] = double(v);
	}

	Vector serial_res(bf.NumRows());
	bf.Mult(serial_gf, serial_res);

	auto parbf = ParBilinearForm(&pfes);
	parbf.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
	parbf.AddInteriorFaceIntegrator(new maxwell::mfemExtension::MaxwellDGOneNormalJumpIntegrator({maxwell::X}, 1.0));
	parbf.Assemble();
	parbf.Finalize();

	ParGridFunction pgf(&pmesh, &serial_gf, partitioning);
	pgf.ExchangeFaceNbrData();
	auto f_nbr_gf = pgf.FaceNbrData();

	Vector local_and_nbr(pgf.Size() + f_nbr_gf.Size());
	for (auto v = 0; v < pgf.Size(); v++){
		local_and_nbr[v] = pgf[v];
	}
	for (auto v = 0; v < f_nbr_gf.Size(); v++){
		local_and_nbr[v+pgf.Size()] = f_nbr_gf[v];
	}

	Vector res_par(parbf.NumRows());
	parbf.Mult(local_and_nbr, res_par);

	Vector serial_cut;
	serial_cut.SetSize(global_dof_indices.size());
	for (int i = 0; i < global_dof_indices.size(); i++) {
    	serial_cut[i] = serial_res[global_dof_indices[i]];
	}

	assert(res_par.Size() == serial_cut.Size());
	for (auto i = 0; i < res_par.Size(); i++){
		EXPECT_EQ(res_par.GetData()[i], serial_cut.GetData()[i]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

TEST_F(MPIFiniteElementSpaceTest, CustomPartitioning) 
{
	auto mesh = Mesh::MakeCartesian1D(10);
	auto partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());

	std::vector<int> values(mesh.GetNE());
	for (auto v = 0; v < mesh.GetNE(); v++){
		values[v] = partitioning[v];
	}

}

}
