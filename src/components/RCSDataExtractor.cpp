#include "RCSDataExtractor.h"

namespace maxwell {

RCSDataExtractor::RCSDataExtractor(const std::string data_folder, const std::string case_name)
{
	// We parse the json file to get info about the simulation data.
	auto case_data = driver::parseJSONfile(maxwellCase(case_name));

	// We assemble the full FES of the original simulation. This is done to ensure the mesh is the original one that 
	// was used during the simulation, to ensure high order meshes if it were necessary.
	auto full_mesh = Mesh::LoadFromFile(driver::assembleMeshString(case_data["model"]["filename"]), 1, 0);
	auto fec = new DG_FECollection(case_data["solver_options"]["order"], full_mesh.Dimension(), BasisType::GaussLobatto);
	auto full_fes = FiniteElementSpace(&full_mesh, fec);

	// This prepares a NTFF submesher that will cut the 'crown mesh' defined by the elements immediately
	// adjacent to the NTFF surface in the scattered field region.
	Probes probes = driver::buildProbes(case_data);
	NearToFarFieldSubMesher ntff_sub(full_mesh, full_fes, buildSurfaceMarker(probes.nearFieldProbes[0].tags, full_fes));
	Mesh crown_mesh(static_cast<Mesh>(*ntff_sub.getSubMesh()));
	auto crown_fes = FiniteElementSpace(&crown_mesh, fec);

	// In order to extract x, y and z node information from MFEM meshes/FES, we need to create a temporary vdim 3 FES
	// and assign a GF to it, then use the crown mesh to get the node coordinate information from them.
	// The nodes are ordered such as X0, X1, X2... Y0, Y1, Y2... Z0, Z1, Z2...
	auto crown_fes_vdim3 = FiniteElementSpace(&crown_mesh, fec, 3);
	GridFunction nodepos(&crown_fes_vdim3);
	crown_mesh.GetNodes(nodepos);
	auto nodepos_dimension_size = nodepos.Size() / 3;
	
	// We initialise our storage vectors to load the values as we iterate through the faces.
	std::vector<double> x, y, z, vx, vy, vz, nx, ny, nz, Ex, Ey, Ez, Hx, Hy, Hz;
	for (int f = 0; f < crown_fes.GetNBE(); f++) {
		//Only the NTFF surface faces are marked with this specific boundary attribute.
		if (crown_fes.GetBdrAttribute(f) == (static_cast<int>(BdrCond::NearToFarField) - 1))
		{

			// This takes care of positions of the nodes.
			Array<int> facedofs;
			crown_fes.GetBdrElementDofs(f, facedofs);
			for (int d = 0; d < facedofs.Size(); d++) {
				x.push_back(nodepos[facedofs[f]]);
				y.push_back(nodepos[facedofs[f] + nodepos_dimension_size]);
				z.push_back(nodepos[facedofs[f] + nodepos_dimension_size * 2]);
			}

			// This takes care of the positions of the geometrical vertices.
			auto be_trans{ crown_mesh.GetBdrElementTransformation(f) };
			auto face{ crown_mesh.GetFace(crown_mesh.GetBdrElementFaceIndex(f)) };
			auto v0{ crown_mesh.GetVertex(face->GetVertices()[0]) };
			auto v1{ crown_mesh.GetVertex(face->GetVertices()[1]) };
			auto v2{ crown_mesh.GetVertex(face->GetVertices()[2]) };
			vx.push_back(v0[0]);
			vx.push_back(v1[0]);
			vx.push_back(v2[0]);
			vy.push_back(v0[1]);
			vy.push_back(v1[1]);
			vy.push_back(v2[1]);
			vz.push_back(v0[2]);
			vz.push_back(v1[2]);
			vz.push_back(v2[2]);

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

			// This takes care of reading the exported GridFunction files at each specific time and
			// then add the field value at the specific degrees of freedom regarding the facedofs
			// to the field vector it belongs to in the same order we organise our facedofs in the previous steps.
			std::vector<std::string> fields({ "/Ex.gf", "/Ey.gf", "/Ez.gf", "/Hx.gf", "/Hy.gf", "/Hz.gf" });
			for (const auto& field : fields) {
				GridFunction A;
				for (auto const& dir_entry : std::filesystem::directory_iterator(data_folder)) {
					if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh" &&
						dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 3) != "rcs" &&
						dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 8) != "farfield") {
						A = getGridFunction(crown_mesh, dir_entry.path().generic_string() + field);
					}
					if (field == "/Ex.gf") {
						for (int d = 0; d < facedofs.Size(); d++) {
							Ex.push_back(A[d]);
						}
					}
					if (field == "/Ey.gf") {
						for (int d = 0; d < facedofs.Size(); d++) {
							Ey.push_back(A[d]);
						}
					}
					if (field == "/Ez.gf") {
						for (int d = 0; d < facedofs.Size(); d++) {
							Ez.push_back(A[d]);
						}
					}
					if (field == "/Hx.gf") {
						for (int d = 0; d < facedofs.Size(); d++) {
							Hx.push_back(A[d]);
						}
					}
					if (field == "/Hy.gf") {
						for (int d = 0; d < facedofs.Size(); d++) {
							Hy.push_back(A[d]);
						}
					}
					if (field == "/Hz.gf") {
						for (int d = 0; d < facedofs.Size(); d++) {
							Hz.push_back(A[d]);
						}
					}
				}
			}

		}
	}
	
	// This takes care of assembling the time vector, by reading through the time.txt files in the data folders.
	std::vector<double> time{ buildTimeVector(data_folder) };

	// This takes care of assembling the incoming field values for the incident planewave.
	auto planewave_data{ buildPlaneWaveData(case_data) };
	Gaussian gauss{ planewave_data.mean, mfem::Vector({-planewave_data.delay}) };
	auto sources = driver::buildSources(case_data);
	std::vector<double> ExInc(time.size()) , EyInc(time.size()), EzInc(time.size());
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
	std::string path(data_folder + "/RCS_LUIS.dat");
	myfile.open(path);
	if (myfile.is_open()) {
		myfile << "N:\n";
		myfile << std::to_string(fec->GetOrder()) + "\n";
		myfile << "NODETOL:\n";
		myfile << std::to_string(1e-10) + "\n";
		myfile << "Nm\n";
		assert(nx.size() == ny.size() == nz.size());
		myfile << "nx ny nz:\n";
		for (int n = 0; n < nx.size(); n++) {
			myfile << std::to_string(nx[n]) + " " + std::to_string(ny[n]) + " " + std::to_string(nz[n]) + "\n";
		}
		assert(x.size() == y.size() == z.size());
		myfile << "x y z:\n";
		for (int n = 0; n < x.size(); n++) {
			myfile << std::to_string(x[n]) + " " + std::to_string(y[n]) + " " + std::to_string(z[n]) + "\n";
		}
		assert(vx.size() == vy.size() == vz.size());
		myfile << "vx vy vz:\n";
		for (int n = 0; n < vx.size(); n++) {
			myfile << std::to_string(vx[n]) + " " + std::to_string(vy[n]) + " " + std::to_string(vz[n]) + "\n";
		}
		myfile << "END HEADER\n";

		myfile << "RCSSTEP:\n";
		myfile << std::to_string(time.size());
		myfile << "END RCSSTEP\n";

		for (int t = 0; t < time.size(); t++) { // We need to correct the time to SI with c.
			myfile << "time:\n";
			myfile << std::to_string(time[t] / physicalConstants::speedOfLight_SI) + "\n";

			assert(ExInc.size() == EyInc.size() == EzInc.size());
			myfile << "ExInc EyInc EzInc:\n";
			myfile << std::to_string(ExInc[t]) + " " + std::to_string(EyInc[t]) + " " + std::to_string(EzInc[t]) + "\n";


			assert(Ex.size() == Ey.size() == Ez.size());
			myfile << "Ex Ey Ez:\n";
			myfile << std::to_string(Ex[t]) + " " + std::to_string(Ey[t]) + " " + std::to_string(Ez[t]) + "\n";


			assert(Hx.size() == Hy.size() == Hz.size());
			myfile << "Hx Hy Hz:\n";
			myfile << std::to_string(Hx[t] / physicalConstants::freeSpaceImpedance_SI) + " " + std::to_string(Hy[t] / physicalConstants::freeSpaceImpedance_SI) + " " + std::to_string(Hz[t] / physicalConstants::freeSpaceImpedance_SI) + "\n";
		}

		myfile.close();

		// FIN
	}
	else {
		throw std::runtime_error("Could not open file to write RCS_LUIS data file.");
	}


}

}