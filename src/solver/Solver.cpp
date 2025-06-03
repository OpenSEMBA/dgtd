#include "Solver.h"
#include <filesystem>
#include <fstream>

using namespace mfem;

namespace maxwell {

std::unique_ptr<FiniteElementSpace> buildFiniteElementSpace(Mesh* m, FiniteElementCollection* fec)
{
#ifdef SEMBA_DGTD_ENABLE_MPI	
	if (dynamic_cast<ParMesh*>(m) != nullptr) {
		auto pm{ dynamic_cast<ParMesh*>(m) };
		return std::make_unique<ParFiniteElementSpace>(pm, fec);
	}
#endif
	if (dynamic_cast<Mesh*>(m) != nullptr) {
		return std::make_unique<FiniteElementSpace>(m, fec);
	}
	throw std::runtime_error("Invalid mesh to build FiniteElementSpace");
}

std::unique_ptr<TimeDependentOperator> Solver::assignEvolutionOperator()
{
	if (opts_.evolution.op == EvolutionOperatorType::Hesthaven) {
		return std::make_unique<HesthavenEvolution>(*fes_, model_, sourcesManager_, opts_.evolution);
	}
	else if (opts_.evolution.op == EvolutionOperatorType::Global) {
		return std::make_unique<GlobalEvolution>(*fes_, model_, sourcesManager_, opts_.evolution);
	}
	else if (opts_.evolution.op == EvolutionOperatorType::Maxwell){
		ProblemDescription pd(model_, probesManager_.probes, sourcesManager_.sources, opts_.evolution);
		return std::make_unique<MaxwellEvolution>(pd, *fes_, sourcesManager_);
	}
}

Solver::Solver(
	const Model& model,
	const Probes& probes,
	const Sources& sources,
	const SolverOptions& options) :
	opts_{ options },
	model_{ model },
	fec_{ opts_.evolution.order, model_.getMesh().Dimension(), opts_.basisType},
	fes_{ buildFiniteElementSpace(& model_.getMesh(), &fec_) },
	fields_{ *fes_ },
	sourcesManager_{ sources, *fes_, fields_ },
	probesManager_ { probes , *fes_, fields_, opts_ },
	time_{0.0}
{

#ifdef ENABLE_STATISTICS_RECORD
	auto initStartTime = std::chrono::steady_clock::now();
	if (!std::filesystem::is_directory("StatisticsExport/") || !std::filesystem::exists("StatisticsExport/")) {
		std::filesystem::create_directory("StatisticsExport/");
	}
	std::filesystem::path statspath("StatisticsExport/" + model_.meshName_ + "/");
	if (std::filesystem::is_directory(statspath) || std::filesystem::exists(statspath)) {
		std::filesystem::remove_all(statspath);
	}
	std::filesystem::create_directory(statspath);
	
#endif

	checkOptionsAreValid(opts_);

	if (opts_.evolution.spectral == true) {
		performSpectralAnalysis(*fes_.get(), model_, opts_.evolution);
	}

	evolTDO_ = assignEvolutionOperator();
	evolTDO_->SetTime(time_);

	if (opts_.timeStep == 0.0) {
		dt_ = estimateTimeStep(model_, opts_, *fes_, evolTDO_.get());
	}
	else {
		dt_ = opts_.timeStep;
	}

	odeSolver_->Init(*evolTDO_);

	probesManager_.updateProbes(time_);


#ifdef ENABLE_STATISTICS_RECORD
	auto initEndTime = std::chrono::steady_clock::now();
	std::ofstream myfile;
	std::string path("StatisticsExport/" + model_.meshName_ + "/" + "statistics.dat");
	myfile.open(path, std::ios::app);
	if (myfile.is_open()) {
		auto runtime = std::chrono::duration<double>(initEndTime - initStartTime).count();
		myfile << std::scientific << std::setprecision(5);
		myfile << "Initialization Time: " << runtime << " (s)\n";
		myfile << std::defaultfloat;
	}
	myfile.close();
#endif

}

void Solver::checkOptionsAreValid(const SolverOptions& opts) const
{
	
	if ((opts.evolution.order < 0) ||
		(opts.finalTime < 0)) {
		throw std::runtime_error("Incorrect parameters in Options");
	}

	if (opts.cfl <= 0.0) {
		throw std::runtime_error("CFL must be positive");
	}

}

const FieldProbe& Solver::getPointProbe(const std::size_t probe) const 
{ 
	return probesManager_.getPointProbe(probe); 
}

const PointProbe& Solver::getFieldProbe(const std::size_t probe) const
{
	return probesManager_.getFieldProbe(probe);
}

double getMinimumInterNodeDistance(FiniteElementSpace& fes)
{
	GridFunction nodes(&fes);
	fes.GetMesh()->GetNodes(nodes);
	double res{ std::numeric_limits<double>::max() };
	for (int e = 0; e < fes.GetMesh()->ElementToElementTable().Size(); ++e) {
		Array<int> dofs;
		fes.GetElementDofs(e, dofs);
		if (dofs.Size() == 1) {
			res = std::min(res, fes.GetMesh()->GetElementSize(e));
		}
		else {
			for (int i = 0; i < dofs.Size(); ++i) {
				for (int j = i + 1; j < dofs.Size(); ++j) {
					res = std::min(res, std::abs(nodes[dofs[i]] - nodes[dofs[j]]));
				}
			}
		}
	}
	return res;
}

bool checkIfElemTypeInMesh(const Mesh& mesh, const Element::Type& type)
{
	for (int e = 0; e < mesh.GetNE(); ++e) {
		if (mesh.GetElementType(e) == type) {
			return true;
		}
	}
	return false;
}


Vector getTimeStepScale(Mesh& mesh)
{
	Vector vol(mesh.GetNE()), dtscale(mesh.GetNE());
	for (int e = 0; e < mesh.GetNE(); ++e) {
		auto el{ mesh.GetElement(e) };
		Vector areasum(mesh.GetNumFaces());
		areasum = 0.0;
		for (int f = 0; f < mesh.GetElement(e)->GetNEdges(); ++f) {
			ElementTransformation* T{ mesh.GetFaceTransformation(f)};
			const IntegrationRule& ir = IntRules.Get(T->GetGeometryType(), T->OrderJ());
			for (int p = 0; p < ir.GetNPoints(); p++)
			{
				const IntegrationPoint& ip = ir.IntPoint(p);
				areasum(e) += ip.weight * T->Weight();
			}
		}
		vol(e) = mesh.GetElementVolume(e);
		dtscale(e) = vol(e) / (areasum(e) / 2.0);
	}
	return dtscale;
}

double getJacobiGQ_RMin(const int order) {
	auto mesh{ Mesh::MakeCartesian1D(1, 2.0) };
	DG_FECollection fec{ order,1,BasisType::GaussLobatto };
	FiniteElementSpace fes{ &mesh, &fec };

	GridFunction nodes(&fes);
	mesh.GetNodes(nodes);

	return std::abs(nodes(0) - nodes(1));
}
std::vector<Source::Position> getVerticesCoordsForElem(const FiniteElementSpace& fes, const ElementId& e, const std::vector<Source::Position>& positions)
{
	Array<int> vertices;
	fes.GetElementVertices(e, vertices);
	std::vector<Source::Position> res(vertices.Size());
	for (auto v{ 0 }; v < vertices.Size(); v++) {
		res[v] = positions[vertices[v]];
	}
	return res;
}

double getSideLength(const Source::Position& va, const Source::Position& vb) {
	return std::sqrt((vb[0] - va[0]) * (vb[0] - va[0]) + (vb[1] - va[1]) * (vb[1] - va[1]));
}

double getElementPerimeter(const std::vector<Source::Position>& vertCoords)
{
	double res = 0.0;
	int n = vertCoords.size();
	for (int i = 0; i < n; i++) {
		res += getSideLength(vertCoords[i], vertCoords[(i + 1) % n]);
	}
	return res;
}

double calc2DMeshTimeStep(FiniteElementSpace& fes)
{
	Vector dtscale{ getTimeStepScale(*fes.GetMesh()) };
	double rmin{ getJacobiGQ_RMin(fes.FEColl()->GetOrder()) };
	auto dt{ dtscale.Min() * rmin * 2.0 / 3.0 / physicalConstants::speedOfLight };
	dt *= 0.75; // Purely heuristic.
	if (checkIfElemTypeInMesh(*fes.GetMesh(),Element::Type::QUADRILATERAL)) {
		return dt / 2.0; // This is purely heuristic.
	}
	else {
		return dt;
	}
}

double calc3DMeshTimeStep(FiniteElementSpace& fes)
{
	Vector dtscale{ getTimeStepScale(*fes.GetMesh()) };
	double rmin{ getJacobiGQ_RMin(fes.FEColl()->GetOrder()) };
	auto dt{ dtscale.Min() * rmin * 2.0 / 3.0 / physicalConstants::speedOfLight };
	dt *= 0.75; // Purely heuristic.
	if (checkIfElemTypeInMesh(*fes.GetMesh(), Element::Type::HEXAHEDRON)) {
		return dt / 6.0; // This is purely heuristic.
	}
	else {
		return dt;
	}
}

double estimateTimeStep(const Model& model, const SolverOptions& opts, FiniteElementSpace& fes, const TimeDependentOperator* tdo)
{
	if (opts.evolution.op != EvolutionOperatorType::Hesthaven) {
		if (model.getConstMesh().Dimension() == 1) {
			double maxTimeStep{ 0.0 };
			if (opts.evolution.order == 0) {
				maxTimeStep = getMinimumInterNodeDistance(fes) / physicalConstants::speedOfLight;
			}
			else {
				maxTimeStep = getMinimumInterNodeDistance(fes) / pow(opts.evolution.order, 1.5) / physicalConstants::speedOfLight;
			}
			return opts.cfl * maxTimeStep;
		}
		else if (model.getConstMesh().Dimension() == 2) {
			return calc2DMeshTimeStep(fes) * opts.cfl;
		}
		else if (model.getConstMesh().Dimension() == 3) {
			return calc3DMeshTimeStep(fes) * opts.cfl / 0.8; // 0.8 is purely heuristic, adjusted from Hesthaven ATS value.
		}
		else {
			throw std::runtime_error("Automatic Time Step Estimation not available for the set dimension.");
		}
	}
	else {
		if (model.getConstMesh().Dimension() == 2) {
			return calc2DMeshTimeStep(fes) * opts.cfl;
		}
		else if (model.getConstMesh().Dimension() == 3) {
			auto maxFscaleVal{ 0.0 };
			const auto& evol = dynamic_cast<const HesthavenEvolution*>(tdo);
			for (auto e{ 0 }; e < model.getConstMesh().GetNE(); e++) {
				const auto& fscaleMax{ evol->getHesthavenElement(e).fscale.maxCoeff() };
				if (fscaleMax > maxFscaleVal) {
					maxFscaleVal = fscaleMax;
				}
			}
			const auto& order = fes.FEColl()->GetOrder();
			return 1.0 * opts.cfl / (maxFscaleVal * order * order);
		}
		else {
			throw std::runtime_error("Automatic Time Step Estimation not available for the set dimension.");
		}
	}
}

#ifdef SHOW_TIMER_INFORMATION
void printSimulationInformation(const double time, const double dt, const double finalTime)
{
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Information is updated every 30 seconds." << std::endl;
	std::cout << "Current Step: " + std::to_string(int(time / dt)) << std::endl;
	std::cout << "Steps Left  : " + std::to_string(int((finalTime - time) / dt)) << std::endl;
	std::cout << std::endl;
	std::cout << "Final Time  : " + std::to_string(finalTime / physicalConstants::speedOfLight_SI * 1e9) + " ns." << std::endl;
	std::cout << "Current Time: " + std::to_string(time / physicalConstants::speedOfLight_SI * 1e9) + " ns." << std::endl;
	std::cout << "Time Step   : " + std::to_string(dt / physicalConstants::speedOfLight_SI * 1e9) + " ns." << std::endl;
	std::cout << std::endl;
}
#endif

void Solver::run()
{

#ifdef ENABLE_STATISTICS_RECORD
	auto runStartTime = std::chrono::steady_clock::now();
#endif

#ifdef SHOW_TIMER_INFORMATION
	auto lastPrintTime{ std::chrono::steady_clock::now() };
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << "-------------SOLVER RUN INFORMATION-------------" << std::endl;
	printSimulationInformation(time_, dt_, opts_.finalTime);
#endif

	while (time_ <= opts_.finalTime - 1e-8*dt_) {

		step();

#ifdef SHOW_TIMER_INFORMATION
		auto currentTime = std::chrono::steady_clock::now();
		if (std::chrono::duration_cast<std::chrono::seconds>
			(currentTime - lastPrintTime).count() >= 30.0) 
		{
			printSimulationInformation(time_, dt_, opts_.finalTime);
			lastPrintTime = currentTime;
		}
		if (this->fields_.getNorml2() > 1e20){
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "Simulation is potentially unstable, verify manually and lower time step." << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
		}

#endif
	}

#ifdef ENABLE_STATISTICS_RECORD
	auto runEndTime = std::chrono::steady_clock::now();
	std::ofstream myfile;
	std::string path("StatisticsExport/" + model_.meshName_ + "/" + "statistics.dat");
	myfile.open(path, std::ios::app);
	if (myfile.is_open()) {
		auto runtime = std::chrono::duration<double>(runEndTime - runStartTime).count();
		myfile << std::scientific << std::setprecision(5);
		myfile << "Simulation Run Time: " << runtime << " (s)\n";
		myfile << std::defaultfloat;
		myfile << "Final Time: " << std::to_string(opts_.finalTime / physicalConstants::speedOfLight_SI * 1e9) << " (ns)\n";
		myfile << "Time Step: " << std::to_string(dt_ / physicalConstants::speedOfLight_SI * 1e9) << " (ns)\n";
		myfile << "Number of Elements: " << std::to_string(fes_.get()->GetNE()) << "\n";
		myfile << "Number of Degrees of Freedom: " << std::to_string(fes_.get()->GetNDofs()) << "\n";
		if (opts_.evolution.op == Global) {
			auto global = dynamic_cast<GlobalEvolution*>(evolTDO_.get());
			myfile << "Number of Global Total Entries: " << std::to_string(int(std::pow(global->getConstGlobalOperator().Size(), 2.0))) << "\n";
			myfile << "Number of Global Non-Zero Entries: " << std::to_string(global->getConstGlobalOperator().NumNonZeroElems()) << "\n";
		}
	}
	myfile.close();
#endif
}

void Solver::step()
{
	double truedt{ std::min(dt_, opts_.finalTime - time_) };
	odeSolver_->Step(fields_.allDOFs(), time_, truedt);
	probesManager_.updateProbes(time_);
}


GeomTagToBoundary Solver::assignAttToBdrByDimForSpectral(Mesh& submesh)
{
	switch (submesh.Dimension()) {
	case 1:
		return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA} };
	case 2:
		switch (submesh.GetElementType(0)) {
		case Element::TRIANGLE:
			return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA } };
		case Element::QUADRILATERAL:
			return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA} };
		default:
			throw std::runtime_error("Incorrect element type for 2D spectral AttToBdr assignation.");
		}
	case 3:
		switch (submesh.GetElementType(0)) {
		case Element::TETRAHEDRON:
			return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA} };
		case Element::HEXAHEDRON:
			return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA}, {5, BdrCond::SMA }, {6, BdrCond::SMA} };
		default:
			throw std::runtime_error("Incorrect element type for 3D spectral AttToBdr assignation.");
		}
	default:
		throw std::runtime_error("Dimension is incorrect for spectral AttToBdr assignation.");
	}

}

Eigen::SparseMatrix<double> Solver::assembleSubmeshedSpectralOperatorMatrix(Mesh& submesh, const FiniteElementCollection& fec, const EvolutionOptions& opts)
{
	Model submodel(submesh, GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo(assignAttToBdrByDimForSpectral(submesh), GeomTagToInteriorBoundary{}));
	FiniteElementSpace subfes(&submesh, &fec);
	Eigen::SparseMatrix<double> local;
	auto numberOfFieldComponents = 2;
	auto numberofMaxDimensions = 3;
	local.resize(numberOfFieldComponents * numberofMaxDimensions * subfes.GetNDofs(), 
		numberOfFieldComponents * numberofMaxDimensions * subfes.GetNDofs());

	EvolutionOptions localopts(opts);
	ProblemDescription pd(submodel, probesManager_.probes, sourcesManager_.sources, localopts);
	DGOperatorFactory dgops(pd, subfes);
	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildDerivativeSubOperator(y), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,z }, -1.0); // MS
		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildDerivativeSubOperator(z), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,y });
		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildDerivativeSubOperator(y), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,z });
		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildDerivativeSubOperator(z), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,y }, -1.0);

		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildOneNormalSubOperator(E, { y }), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,z }); // MFN
		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildOneNormalSubOperator(E, { z }), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,y }, -1.0);
		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildOneNormalSubOperator(H, { y }), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,z }, -1.0);
		allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildOneNormalSubOperator(H, { z }), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,y });

		if (opts.alpha > 0.0) {

			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildZeroNormalSubOperator(H), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { x }, -1.0); // MP
			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildZeroNormalSubOperator(E), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { x }, -1.0);

			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildTwoNormalSubOperator(H, { X, x }), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { X,x }); //MPNN
			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildTwoNormalSubOperator(H, { Y, x }), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { Y,x });
			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(H), *dgops.buildTwoNormalSubOperator(H, { Z, x }), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { Z,x });
			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildTwoNormalSubOperator(E, { X, x }), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { X,x });
			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildTwoNormalSubOperator(E, { Y, x }), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { Y,x });
			allocateDenseInEigen(buildByMult(*dgops.buildInverseMassMatrixSubOperator(E), *dgops.buildTwoNormalSubOperator(E, { Z, x }), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { Z,x });

		}

	}
	return local;
}

double Solver::findMaxEigenvalueModulus(const Eigen::VectorXcd& eigvals)
{
	auto res{ 0.0 };
	for (int i = 0; i < eigvals.size(); ++i) {
		auto modulus{ sqrt(pow(eigvals[i].real(),2.0) + pow(eigvals[i].imag(),2.0)) };
		if (modulus <= 1.0 && modulus >= res) {
			res = modulus;
		}
	}
	return res;
}

void reassembleSpectralBdrForSubmesh(SubMesh* submesh) 
{
	switch (submesh->GetElementType(0)) {
	case Element::SEGMENT:
		for (int i = 0; i < submesh->GetParentVertexIDMap().Size(); ++i) {
			submesh->AddBdrPoint(i, i + 1);
		}
		submesh->FinalizeMesh();
		break;
	case Element::TRIANGLE:
		for (int i = 0; i < submesh->GetNBE(); ++i) {
			submesh->SetBdrAttribute(i, i + 1);
		}
		submesh->FinalizeMesh();
		break;
	case Element::QUADRILATERAL:
		for (int i = 0; i < submesh->GetNBE(); ++i) {
			submesh->SetBdrAttribute(i, i + 1);
		}
		submesh->FinalizeMesh();
		break;
	case Element::TETRAHEDRON:
		for (int i = 0; i < submesh->GetNBE(); ++i) {
			submesh->SetBdrAttribute(i, i + 1);
		}
		submesh->FinalizeMesh();
		break;
	case Element::HEXAHEDRON:
		for (int i = 0; i < submesh->GetNBE(); ++i) {
			submesh->SetBdrAttribute(i, i + 1);
		}
		submesh->FinalizeMesh();
		break;
	default:
		throw std::runtime_error("Incorrect element type for Bdr Spectral assignation.");
	}
}

void Solver::evaluateStabilityByEigenvalueEvolutionFunction(
	Eigen::VectorXcd& eigenvals, 
	MaxwellEvolution& maxwellEvol)
{
	auto real { toMFEMVector(eigenvals.real()) };
	auto realPre = real;
	auto imag { toMFEMVector(eigenvals.imag()) };
	auto imagPre = imag;
	auto time { 0.0 };
	maxwellEvol.SetTime(time);
	odeSolver_->Init(maxwellEvol);
	odeSolver_->Step(real, time, opts_.timeStep);
	time = 0.0;
	maxwellEvol.SetTime(time);
	odeSolver_->Init(maxwellEvol);
	odeSolver_->Step(imag, time, opts_.timeStep);
	
	for (int i = 0; i < real.Size(); ++i) {
		
		auto modPre{ sqrt(pow(realPre[i],2.0) + pow(imagPre[i],2.0)) };
		auto mod   { sqrt(pow(real[i]   ,2.0) + pow(imag[i]   ,2.0)) };

		if (modPre != 0.0) {
			if (mod / modPre > 1.0) {
				throw std::runtime_error("The coefficient between the modulus of a time evolved eigenvalue and its original value is higher than 1.0 - RK4 instability.");
			}
		}
	}
}

void Solver::performSpectralAnalysis(const FiniteElementSpace& fes, Model& model, const EvolutionOptions& opts)
{
	Array<int> domainAtts(1);
	domainAtts[0] = 501;
	auto mesh{ model.getConstMesh() };
	auto meshCopy{ mesh };

	for (int elem = 0; elem < meshCopy.GetNE(); ++elem) {

		auto preAtt(meshCopy.GetAttribute(elem));
		meshCopy.SetAttribute(elem, domainAtts[0]);
		auto submesh{ SubMesh::CreateFromDomain(meshCopy,domainAtts) };
		meshCopy.SetAttribute(elem, preAtt);
		submesh.SetAttribute(0, preAtt);

		reassembleSpectralBdrForSubmesh(&submesh);

		auto eigenvals{ 
			assembleSubmeshedSpectralOperatorMatrix(submesh, *fes.FEColl(), opts).toDense().eigenvalues() 
		};
		FiniteElementSpace submeshFES{ &submesh, fes.FEColl() };
		Model model{ submesh,
			GeomTagToMaterialInfo{},
			GeomTagToBoundaryInfo(assignAttToBdrByDimForSpectral(submesh),GeomTagToInteriorBoundary{})
		};
		SourcesManager srcs{ Sources(), submeshFES, fields_ };
		ProblemDescription pd(model, probesManager_.probes, sourcesManager_.sources, opts_.evolution);
		MaxwellEvolution evol(pd, submeshFES, sourcesManager_);
		evaluateStabilityByEigenvalueEvolutionFunction(eigenvals, evol);
	}
}


}
