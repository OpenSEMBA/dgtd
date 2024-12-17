#include "Evolution.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

void evalConductivity(const Vector& cond, const Vector& in, Vector& out)
{
	for (auto v{ 0 }; v < cond.Size(); v++) {
		out[v] -= cond[v] * in[v];
	}
}

void changeSignOfFieldGridFuncs(FieldGridFuncs& gfs)
{
	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			gfs[f][d] *= -1.0;
		}
	}
}

const FieldGridFuncs evalTimeVarFunction(const Time time, SourcesManager& sm)
{
	auto res{ sm.evalTimeVarField(time, sm.getGlobalTFSFSpace()) };
	auto func_g_sf = res;
	sm.markDoFSforTForSF(res, true);
	{
		if (sm.getTFSFSubMesher().getSFSubMesh() != NULL) {
			sm.markDoFSforTForSF(func_g_sf, false);
			for (int f : {E, H}) {
				for (int x{ 0 }; x <= Z; x++) {
					res[f][x] -= func_g_sf[f][x];
					res[f][x] *= 0.5;
				}
			}
		}
	}
	return res;
}

Evolution::Evolution(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{
#ifdef SHOW_TIMER_INFORMATION
	auto startTime{ std::chrono::high_resolution_clock::now() };
#endif


#ifdef SHOW_TIMER_INFORMATION
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << std::endl;
#endif

	if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model.getTotalFieldScatteredFieldToMarker().end()) {
		srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));
		auto globalTFSFfes{ srcmngr_.getGlobalTFSFSpace() };
		Model modelGlobal = Model(*globalTFSFfes->GetMesh(), GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Assembling TFSF Inverse Mass Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			MInvTFSF_[f] = buildInverseMassMatrix(f, modelGlobal, *globalTFSFfes);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono:: >
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass Zero-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			MP_GTFSF_[f] = buildByMult(*MInvTFSF_[f], *buildZeroNormalOperator(f, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass One-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					MFN_GTFSF_[f][f2][d] = buildByMult(*MInvTFSF_[f], *buildOneNormalOperator(f2, {d}, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
				}
			}
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass Two-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						MFNN_GTFSF_[f][f2][d][d2] = buildByMult(*MInvTFSF_[f], *buildTwoNormalOperator(f2, {d, d2}, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
					}
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		MInv_[f] = buildInverseMassMatrix(f, model_, fes_);
	}

	if (model_.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			MPB_[f] = buildByMult(*MInv_[f], *buildZeroNormalIBFIOperator(f, model_, fes_, opts_), fes_);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					MFNB_[f][f2][d] = buildByMult(*MInv_[f], *buildOneNormalIBFIOperator(f2, { d }, model_, fes_, opts_), fes_);
				} 
			}
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						MFNNB_[f][f2][d][d2] = buildByMult(*MInv_[f], *buildTwoNormalIBFIOperator(f2, { d, d2 }, model_, fes_, opts_), fes_);
					}
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			MS_[f][d] = buildByMult(*MInv_[f], *buildDerivativeOperator(d, fes_), fes_);
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		MP_[f] = buildByMult(*MInv_[f], *buildZeroNormalOperator(f, model_, fes_, opts_), fes_);
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			for (auto f2 : { E, H }) {
				MFN_[f][f2][d] = buildByMult(*MInv_[f], *buildOneNormalOperator(f2, { d }, model_, fes_, opts_), fes_);
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			for (auto f2 : { E, H }) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					MFNN_[f][f2][d][d2] = buildByMult(*MInv_[f], *buildTwoNormalOperator(f2, { d, d2 }, model_, fes_, opts_), fes_);
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Operator assembly finished" << std::endl;
	std::cout << std::endl;
#endif

	CND_ = buildConductivityCoefficients(model_, fes_);
 }

void Evolution::Mult(const Vector& in, Vector& out) const
{
	const auto& dim{ fes_.GetMesh()->Dimension() };

	std::array<Vector, 3> eOld, hOld;
	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetDataAndSize(in.GetData() + d * fes_.GetNDofs(), fes_.GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d + 3) * fes_.GetNDofs(), fes_.GetNDofs());
		eNew[d].SetSpace(&fes_);
		hNew[d].SetSpace(&fes_);
		eNew[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hNew[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
		eNew[d] = 0.0;
		hNew[d] = 0.0;
	}

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;
		
		evalConductivity(CND_, eOld[x], eNew[x]);

		//Centered
		MS_[H][y]		 ->AddMult(eOld[z], hNew[x],-1.0);
		MS_[H][z]		 ->AddMult(eOld[y], hNew[x]);
		MS_[E][y]		 ->AddMult(hOld[z], eNew[x]);
		MS_[E][z]		 ->AddMult(hOld[y], eNew[x],-1.0);
		
		MFN_[H][E][y]  	 ->AddMult(eOld[z], hNew[x], 1.0);
		MFN_[H][E][z]	 ->AddMult(eOld[y], hNew[x],-1.0);
		MFN_[E][H][y]  	 ->AddMult(hOld[z], eNew[x],-1.0);
		MFN_[E][H][z]	 ->AddMult(hOld[y], eNew[x], 1.0);

		if (opts_.fluxType == FluxType::Upwind) {

			MFNN_[H][H][X][x]->AddMult(hOld[X], hNew[x], 1.0);
			MFNN_[H][H][Y][x]->AddMult(hOld[Y], hNew[x], 1.0);
			MFNN_[H][H][Z][x]->AddMult(hOld[Z], hNew[x], 1.0);
			MP_[H]			 ->AddMult(hOld[x], hNew[x],-1.0);
			
			MFNN_[E][E][Y][x]->AddMult(eOld[Y], eNew[x], 1.0);
			MFNN_[E][E][X][x]->AddMult(eOld[X], eNew[x], 1.0);
			MFNN_[E][E][Z][x]->AddMult(eOld[Z], eNew[x], 1.0);
			MP_[E]			 ->AddMult(eOld[x], eNew[x],-1.0);
		}

		if (model_.getInteriorBoundaryToMarker().size() != 0) {
			
			MFNB_[H][E][y]->AddMult(eOld[z], hNew[x]);
			MFNB_[H][E][z]->AddMult(eOld[y], hNew[x], -1.0);
			MFNB_[E][H][y]->AddMult(hOld[z], eNew[x], -1.0);
			MFNB_[E][H][z]->AddMult(hOld[y], eNew[x]);

			if (opts_.fluxType == FluxType::Upwind) {

				MFNNB_[H][H][X][x]->AddMult(hOld[X], hNew[x]);
				MFNNB_[H][H][Y][x]->AddMult(hOld[Y], hNew[x]);
				MFNNB_[H][H][Z][x]->AddMult(hOld[Z], hNew[x]);
				MPB_[H]->AddMult(hOld[x], hNew[x], -1.0);

				MFNNB_[E][E][X][x]->AddMult(eOld[X], eNew[x]);
				MFNNB_[E][E][Y][x]->AddMult(eOld[Y], eNew[x]);
				MFNNB_[E][E][Z][x]->AddMult(eOld[Z], eNew[x]);
				MPB_[E]->AddMult(eOld[x], eNew[x], -1.0);
			}
		}
	}

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<Planewave*>(source.get())) {
			
			auto func { evalTimeVarFunction(GetTime(),srcmngr_) };

			std::array<GridFunction, 3> eTemp, hTemp;

			for (int d = X; d <= Z; d++) {
				eTemp[d].SetSpace(srcmngr_.getGlobalTFSFSpace());
				hTemp[d].SetSpace(srcmngr_.getGlobalTFSFSpace());
				eTemp[d] = 0.0;
				hTemp[d] = 0.0;
			}

			for (int x = X; x <= Z; x++) {
				int y = (x + 1) % 3;
				int z = (x + 2) % 3;

				MaxwellTransferMap eMap(eTemp[x], eNew[x]);
				MaxwellTransferMap hMap(hTemp[x], hNew[x]);

				//Centered

				MFN_GTFSF_[H][E][y]->Mult(func[E][z], hTemp[x]);
				eMap.TransferSub(hTemp[x], hNew[x]);
				MFN_GTFSF_[H][E][z]->Mult(func[E][y], hTemp[x]);
				eMap.TransferAdd(hTemp[x], hNew[x]);
				MFN_GTFSF_[E][H][y]->Mult(func[H][z], eTemp[x]);
				eMap.TransferAdd(eTemp[x], eNew[x]);
				MFN_GTFSF_[E][H][z]->Mult(func[H][y], eTemp[x]);
				eMap.TransferSub(eTemp[x], eNew[x]);

				if (opts_.fluxType == FluxType::Upwind) {
					MFNN_GTFSF_[H][H][X][x]->Mult(func[H][X], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MFNN_GTFSF_[H][H][Y][x]->Mult(func[H][Y], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MFNN_GTFSF_[H][H][Z][x]->Mult(func[H][Z], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MP_GTFSF_[H]           ->Mult(func[H][x], hTemp[x]);
					hMap.TransferAdd(hTemp[x], hNew[x]);

					MFNN_GTFSF_[E][E][X][x]->Mult(func[E][X], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MFNN_GTFSF_[E][E][Y][x]->Mult(func[E][Y], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MFNN_GTFSF_[E][E][Z][x]->Mult(func[E][Z], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MP_GTFSF_[E]           ->Mult(func[E][x], eTemp[x]);
					eMap.TransferAdd(eTemp[x], eNew[x]);
				}

			}
		}
	}
}

void HesthavenEvolution::emplaceEmat(const DynamicMatrix& surface_matrix, const FaceId f, HesthavenElement& hestElem)
{
	std::set<DynamicMatrix, MatrixCompareLessThan>::iterator it = matrixStorage_.find(surface_matrix);
	switch (f) {
	case 0:
		hestElem.emat.face0 = &(*it);
		break;
	case 1:
		hestElem.emat.face1 = &(*it);
		break;
	case 2:
		hestElem.emat.face2 = &(*it);
		break;
	default:
		throw std::runtime_error("Hesthaven Evolution only supports triangles for now.");
	}
}

DynamicMatrix assembleInverseMassMatrix(FiniteElementSpace& fes)
{
	BilinearForm bf(&fes);
	ConstantCoefficient one(1.0);
	bf.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
	bf.Assemble();
	bf.Finalize();

	return toEigen(*bf.SpMat().ToDenseMatrix());
}

DynamicMatrix getReferenceInverseMassMatrix(const Mesh& mesh, const int order)
{
	auto m{ Mesh::MakeCartesian2D(1, 1, mesh.GetElementType(0)) };
	auto fec{ L2_FECollection(order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };
	auto mass_mat{ assembleInverseMassMatrix(fes) };

	DynamicMatrix res = getElementMassMatrixFromGlobal(0, mass_mat);
	return res;
}

HesthavenEvolution::HesthavenEvolution(FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& opts) :
	fes_(fes),
	model_(model),
	srcmngr_(srcmngr),
	opts_(opts)
{
	Array<int> elementMarker;
	elementMarker.Append(hesthavenMeshingTag);

	connectivity_ =  assembleGlobalConnectivityMap(model_.getMesh(), dynamic_cast<const L2_FECollection*>(fes.FEColl()));

	for (auto e{ 0 }; e < model.getConstMesh().GetNE(); e++)
	{
		HesthavenElement hestElem;
		hestElem.id = e;
		hestElem.geom = model.getConstMesh().GetElementBaseGeometry(e);

		auto m{ Mesh(model.getMesh()) };

		m.SetAttribute(e, hesthavenMeshingTag);
		auto sm = SubMesh::CreateFromDomain(m, elementMarker);
		auto fec = dynamic_cast<const L2_FECollection*>(fes.FEColl());
		FiniteElementSpace sm_fes(&sm, dynamic_cast<const L2_FECollection*>(fes.FEColl()));

		auto boundary_markers = assembleBoundaryMarkers(sm_fes);
		for (auto f{ 0 }; f < sm_fes.GetNF(); f++) {
			sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
		}

		for (auto f{ 0 }; f < sm.GetNEdges(); f++){
			auto surface_matrix{ assembleConnectivityFaceMassMatrix(sm_fes, boundary_markers[f]) };
			std::set<DynamicMatrix, MatrixCompareLessThan>::iterator it = matrixStorage_.find(surface_matrix);
			if (it == matrixStorage_.end()) {
				matrixStorage_.insert(surface_matrix);
				emplaceEmat(surface_matrix, f, hestElem);
			}
			else {
				emplaceEmat(surface_matrix, f, hestElem);
			}
		}

		for (auto f{ 0 }; f < sm.GetNEdges(); f++) {
			Vector normal(sm.SpaceDimension());
			ElementTransformation* f_trans = sm.GetEdgeTransformation(f);
			f_trans->SetIntPoint(&Geometries.GetCenter(f_trans->GetGeometryType()));
			CalcOrtho(f_trans->Jacobian(), normal);
			hestElem.normals.X.resize(model.getConstMesh().GetElement(e)->GetNEdges() * (sm_fes.FEColl()->GetOrder() + 1));
			hestElem.normals.Y.resize(model.getConstMesh().GetElement(e)->GetNEdges() * (sm_fes.FEColl()->GetOrder() + 1));
			hestElem.fscale   .resize(model.getConstMesh().GetElement(e)->GetNEdges() * (sm_fes.FEColl()->GetOrder() + 1));
			for (auto b{ 0 }; b < sm_fes.FEColl()->GetOrder() + 1; b++) { //hesthaven requires normals to be stored per node at face
				hestElem.normals.X[b] = normal[0];
				hestElem.normals.Y[b] = normal[1];
				hestElem.fscale[b] = abs(normal[0]) + abs(normal[1]); //likewise for fscale, surface per volume ratio per node at face
			}
		}

		hestElemStorage_.push_back(hestElem);

	}


}

void HesthavenEvolution::Mult(const Vector& in, Vector& out)
{
	const Eigen::Map<Eigen::VectorXd> Ex_in(in.GetData() + 0 * fes_.GetNDofs(), fes_.GetNDofs(), 1);
	const Eigen::Map<Eigen::VectorXd> Ey_in(in.GetData() + 1 * fes_.GetNDofs(), fes_.GetNDofs(), 1);
	const Eigen::Map<Eigen::VectorXd> Ez_in(in.GetData() + 2 * fes_.GetNDofs(), fes_.GetNDofs(), 1);
	const Eigen::Map<Eigen::VectorXd> Hx_in(in.GetData() + 3 * fes_.GetNDofs(), fes_.GetNDofs(), 1);
	const Eigen::Map<Eigen::VectorXd> Hy_in(in.GetData() + 4 * fes_.GetNDofs(), fes_.GetNDofs(), 1);
	const Eigen::Map<Eigen::VectorXd> Hz_in(in.GetData() + 5 * fes_.GetNDofs(), fes_.GetNDofs(), 1);

	// ---JUMPS--- //

	Eigen::VectorXd dEx(fes_.GetNDofs()), dEy(fes_.GetNDofs()), dEz(fes_.GetNDofs()),
					dHx(fes_.GetNDofs()), dHy(fes_.GetNDofs()), dHz(fes_.GetNDofs());

	for (auto v{ 0 }; v < connectivity_.size(); v++) {
		dEx[v] = Ex_in[connectivity_[v].first] - Ex_in[connectivity_[v].second];
		dEy[v] = Ey_in[connectivity_[v].first] - Ey_in[connectivity_[v].second];
		dEz[v] = Ez_in[connectivity_[v].first] - Ez_in[connectivity_[v].second];
		dHx[v] = Hx_in[connectivity_[v].first] - Hx_in[connectivity_[v].second];
		dHy[v] = Hy_in[connectivity_[v].first] - Hy_in[connectivity_[v].second];
		dHz[v] = Hz_in[connectivity_[v].first] - Hz_in[connectivity_[v].second];
	}

	// ----------- //

	// Extend to all elements

	auto ndotdH = hestElemStorage_[0].normals.X.asDiagonal() * dHx + hestElemStorage_[0].normals.Y.asDiagonal() * dHy + hestElemStorage_[0].normals.Z.asDiagonal() * dHz;
	auto ndotdE = hestElemStorage_[0].normals.X.asDiagonal() * dEx + hestElemStorage_[0].normals.Y.asDiagonal() * dEy + hestElemStorage_[0].normals.Z.asDiagonal() * dEz;

	auto inverseMassMatrix{ getReferenceInverseMassMatrix(model_.getConstMesh(), fes_.FEColl()->GetOrder()) };



}

}

