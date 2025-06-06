#include "RCSDataExtractor.h"

#include <mfem.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip> 

#include "RCSManager.h"
#include "SubMesher.h"
#include "../test/TestUtils.h"
#include "math/PhysicalConstants.h"
#include "evolution/HesthavenEvolutionMethods.h"

namespace maxwell {

std::unique_ptr<SparseMatrix> assembleNtFFMatrix(FiniteElementSpace& fes) {

	Array<int> n2ffmarker(static_cast<int>(BdrCond::NearToFarField));
	n2ffmarker = 0;
	n2ffmarker[static_cast<int>(BdrCond::NearToFarField) - 1] = 1;
	BilinearForm boundaryForm(&fes);
	boundaryForm.AddBdrFaceIntegrator(new mfemExtension::MaxwellDGZeroNormalJumpIntegrator(1.0), n2ffmarker);
	boundaryForm.Assemble();
	boundaryForm.Finalize();

	return std::make_unique<SparseMatrix>(boundaryForm.SpMat());

}

struct ArrayLess
{
	bool operator()(const Array<int>& a, const Array<int>& b) const
	{
		if (a.Size() != b.Size()) return a.Size() < b.Size();
		for (int i = 0; i < a.Size(); ++i)
		{
			if (a[i] != b[i]) return a[i] < b[i];
		}
		return false;
	}
};

std::map<FaceId, Array<int>> getNodesFromNtFFFaces(FiniteElementSpace& fes) {

	std::map<FaceId, Array<int>> res;

	auto spmat = assembleNtFFMatrix(fes);

	std::set<NodeId> nodes;

	double tol = 1e-10;
	for (auto r = 0; r < spmat->NumRows(); r++) {
		Array<int> cols;
		Vector vals;
		Nodes nonZeroEntries;
		spmat->GetRow(r, cols, vals);
		for (auto v = 0; v < vals.Size(); v++) {
			if (std::abs(vals[v]) > tol) {
				nodes.insert(cols[v]);
			}
		}
	}

	std::vector<Array<int>> tempNodeSet;
	auto it = nodes.begin();
	while (std::distance(it, nodes.end()) >= 3) {
		Array<int> nodesToMap(3);
		for (int i = 0; i < 3; ++i) {
			nodesToMap[i] = *it;
			++it;
		}
		tempNodeSet.push_back(nodesToMap);
	}

	std::vector<FaceId> taggedFaces;
	for (auto b = 0; b < fes.GetNBE(); b++) {
		if (fes.GetBdrAttribute(b) == static_cast<int>(BdrCond::NearToFarField)) {
			taggedFaces.push_back(b);
		}
	}

	for (auto b = 0; b < taggedFaces.size(); b++) {
		res[taggedFaces[b]] = tempNodeSet[b];
	}

	return res;
}

RCSDataExtractor::RCSDataExtractor(const std::string data_folder, const std::string case_name)
{
	// We parse the json file to get info about the simulation data.
	auto case_data = driver::parseJSONfile(maxwellCase(case_name));

	// We assemble the full FES of the original simulation. This is done to ensure the mesh is the original one that 
	// was used during the simulation, to ensure high order meshes if it were necessary.
	auto mesh = Mesh::LoadFromFile(driver::assembleMeshString(case_data["model"]["filename"]), 1, 0);
	auto fec = new DG_FECollection(case_data["solver_options"]["order"], mesh.Dimension(), BasisType::GaussLobatto);

	// This prepares a NTFF submesher that will cut the 'crown mesh' defined by the elements immediately
	// adjacent to the NTFF surface in the scattered field region.
	Probes probes = driver::buildProbes(case_data);
	ParMesh full_mesh = ParMesh(MPI_COMM_WORLD, mesh);
	ParFiniteElementSpace pfull_fes(&full_mesh, fec);
	NearToFarFieldSubMesher ntff_sub(full_mesh, pfull_fes, buildSurfaceMarker(probes.nearFieldProbes[0].tags, pfull_fes));
	ParMesh crown_mesh(static_cast<ParMesh>(*ntff_sub.getSubMesh()));
	auto crown_fes = ParFiniteElementSpace(&crown_mesh, fec);

	// In order to extract x, y and z node information from MFEM meshes/FES, we need to create a temporary vdim 3 FES
	// and assign a GF to it, then use the crown mesh to get the node coordinate information from them.
	// The nodes are ordered such as X0, X1, X2... Y0, Y1, Y2... Z0, Z1, Z2...
	auto crown_fes_vdim3 = ParFiniteElementSpace(&crown_mesh, fec, 3);
	ParGridFunction nodepos(&crown_fes_vdim3);
	crown_mesh.GetNodes(nodepos);
	auto nodepos_dimension_size = nodepos.Size() / 3;

	ParaViewDataCollection pd("crown_mesh", &crown_mesh);
	pd.SetPrefixPath("ParaView");
	const auto order{ 1 };
	pd.SetLevelsOfDetail(1);
	pd.SetHighOrderOutput(false);
	pd.SetDataFormat(VTKFormat::BINARY);
	pd.Save();

	auto nodeIds = getNodesFromNtFFFaces(crown_fes);

	// We initialise our storage vectors to load the values as we iterate through the faces.
	std::vector<double> x, y, z, vx, vy, vz, nx, ny, nz;
	std::vector<std::vector<double>> Ex, Ey, Ez, Hx, Hy, Hz;
	for (const auto& [f, nodes] : nodeIds) {

		// This takes care of positions of the nodes.
		for (int d = 0; d < nodes.Size(); d++) {
			x.push_back(nodepos[nodes[d]]);
			y.push_back(nodepos[nodes[d] + nodepos_dimension_size]);
			z.push_back(nodepos[nodes[d] + nodepos_dimension_size * 2]);
		}

		// This takes care of the positions of the geometrical vertices. This is order 1, so x = vx, y = vy, z = vz.
		auto be_trans{ crown_mesh.GetBdrElementTransformation(f) };
		auto face{ crown_mesh.GetFace(crown_mesh.GetBdrElementFaceIndex(f)) };
		//auto v0{ crown_mesh.GetVertex(face->GetVertices()[0]) }; //These vertices wouldn't be ordered in a similar way to 
		//auto v1{ crown_mesh.GetVertex(face->GetVertices()[1]) }; //the node order
		//auto v2{ crown_mesh.GetVertex(face->GetVertices()[2]) };
		for (int d = 0; d < nodes.Size(); d++) {
			vx.push_back(nodepos[nodes[d]]);
			vy.push_back(nodepos[nodes[d] + nodepos_dimension_size]);
			vz.push_back(nodepos[nodes[d] + nodepos_dimension_size * 2]);
		}

		// This takes care of the normal values, this assumes we're only considering order 1 isoparametric problems AND 
		// that we need to invert the direction of the normal because it points towards the origin by default.
		Vector normals(crown_mesh.Dimension());
		be_trans->SetIntPoint(&Geometries.GetCenter(be_trans->GetGeometryType()));
		CalcOrtho(be_trans->Jacobian(), normals);
		for (auto v = 0; v < face->GetNVertices(); v++) {
			nx.push_back(normals[0] * -1.0);
			ny.push_back(normals[1] * -1.0);
			nz.push_back(normals[2] * -1.0);
		}
	}

	Nodes allNodes;
	for (const auto& [k, nodes] : nodeIds) {
		for (auto node : nodes) {
			allNodes.push_back(node);
		}
	}

	// This takes care of reading the exported GridFunction files at each specific time and
	// then add the field value at the specific degrees of freedom regarding the facedofs
	// to the field vector it belongs to in the same order we organise our facedofs in the previous steps.

	std::vector<std::string> fields({ "/Ex.gf", "/Ey.gf", "/Ez.gf", "/Hx.gf", "/Hy.gf", "/Hz.gf" });

	std::map<std::string, std::vector<ParGridFunction>> GridFuncs;
	for (const auto& field : fields) {
		std::vector<ParGridFunction> A;
		for (auto const& dir_entry : std::filesystem::directory_iterator(data_folder)) {
			if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh" &&
				dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 3) != "rcs" &&
				dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 8) != "farfield") {
				A.push_back(getGridFunction(crown_mesh, dir_entry.path().generic_string() + field));
			}
		}
		GridFuncs[field] = A;
	}

	// This takes care of assembling the time vector, by reading through the time.txt files in the data folders.
	std::vector<double> time{ buildTimeVector(data_folder) };

	for (const auto& field : fields) {
		for (auto g = 0; g < GridFuncs[field].size(); g++) {
			std::vector<double> fieldvals(allNodes.size());
			for (auto v = 0; v < fieldvals.size(); v++) {
				fieldvals[v] = GridFuncs[field][g][allNodes[v]];
			}
			if (field == "/Ex.gf") {
				Ex.push_back(fieldvals);
			}
			if (field == "/Ey.gf") {
				Ey.push_back(fieldvals);
			}
			if (field == "/Ez.gf") {
				Ez.push_back(fieldvals);
			}
			if (field == "/Hx.gf") {
				Hx.push_back(fieldvals);
			}
			if (field == "/Hy.gf") {
				Hy.push_back(fieldvals);
			}
			if (field == "/Hz.gf") {
				Hz.push_back(fieldvals);
			}
		}
	}

	// This takes care of assembling the incoming field values for the incident planewave.
	auto planewave_spread = case_data["sources"][0]["magnitude"]["spread"];
	auto planewave_mean = driver::assemble3DVector(case_data["sources"][0]["magnitude"]["mean"]);
	auto planewave_propagation = driver::assemble3DVector(case_data["sources"][0]["propagation"]);
	auto proj_mean = planewave_mean * planewave_propagation / planewave_propagation.Norml2();
	Gaussian gauss{ planewave_spread, mfem::Vector({proj_mean}) };
	auto sources = driver::buildSources(case_data);
	std::vector<double> ExInc(time.size()), EyInc(time.size()), EzInc(time.size());
	for (const auto& source : sources) {
		if (dynamic_cast<TotalField*>(source.get())) {
			Position pos({ 0.0, 0.0, 0.0 });
			for (int t = 0; t < time.size(); t++) {
				ExInc[t] = source->eval(pos, time[t], E, X);
				EyInc[t] = source->eval(pos, time[t], E, Y);
				EzInc[t] = source->eval(pos, time[t], E, Z);
			}
		}
	}

	// Finally, we write the file following the pattern used by Luis in his code. (dgtd-rcs 2008)


	std::ofstream myfile;
	myfile << std::scientific << std::setprecision(10);
	std::string path("./testData/maxwellInputs/" + case_name + "/RCS_LUIS.dat");
	myfile.open(path);

	if (myfile.is_open()) {
		myfile << "N: ";
		std::ostringstream ss;
		ss << fec->GetOrder();
		myfile << ss.str() + "\n";

		myfile << "NODETOL: ";
		ss.str("");
		ss.clear();
		ss << 1e-10;
		myfile << ss.str() + " Nm \n";

		assert(nx.size() == ny.size() == nz.size());
		myfile << "nx ny nz:\n";
		for (int n = 0; n < nx.size(); n++) {
			ss.str("");
			ss.clear();
			ss << nx[n] << " " << ny[n] << " " << nz[n];
			myfile << ss.str() + "\n";
		}

		assert(x.size() == y.size() == z.size());
		myfile << "x y z:\n";
		for (int n = 0; n < x.size(); n++) {
			ss.str("");
			ss.clear();
			ss << x[n] << " " << y[n] << " " << z[n];
			myfile << ss.str() + "\n";
		}

		assert(vx.size() == vy.size() == vz.size());
		myfile << "vx vy vz:\n";
		for (int n = 0; n < vx.size(); n++) {
			ss.str("");
			ss.clear();
			ss << vx[n] << " " << vy[n] << " " << vz[n];
			myfile << ss.str() + "\n";
		}

		myfile << "END HEADER\n";

		for (int t = 0; t < time.size(); t++) {
			myfile << "RCSSTEP:\n";
			myfile << "time:\n";

			ss.str("");
			ss.clear();
			ss << time[t] / physicalConstants::speedOfLight_SI;
			myfile << ss.str() + "\n";

			assert(ExInc.size() == EyInc.size() == EzInc.size());
			myfile << "ExInc EyInc EzInc:\n";
			ss.str("");
			ss.clear();
			ss << ExInc[t] << " " << EyInc[t] << " " << EzInc[t];
			myfile << ss.str() + "\n";

			assert(Ex.size() == Ey.size() == Ez.size());
			myfile << "Ex Ey Ez:\n";
			for (auto v = 0; v < Ex[t].size(); v++) {
				ss.str("");
				ss.clear();
				ss << Ex[t][v] << " " << Ey[t][v] << " " << Ez[t][v];
				myfile << ss.str() + "\n";
			}

			assert(Hx.size() == Hy.size() == Hz.size());
			myfile << "Hx Hy Hz:\n";
			for (auto v = 0; v < Hx[t].size(); v++) {
				ss.str("");
				ss.clear();
				ss  << Hx[t][v] / physicalConstants::freeSpaceImpedance_SI << " "
					<< Hy[t][v] / physicalConstants::freeSpaceImpedance_SI << " "
					<< Hz[t][v] / physicalConstants::freeSpaceImpedance_SI;
				myfile << ss.str() + "\n";
			}

			myfile << "END RCSSTEP\n";
		}

		myfile.close();

		// FIN
	}
	else {
		throw std::runtime_error("Could not open file to write RCS_LUIS data file.");
	}


}

}