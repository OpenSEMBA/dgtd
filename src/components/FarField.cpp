#include "components/FarField.h"

namespace maxwell {

using namespace mfem;

const Vector buildObsPointVec(const SphericalAngles& angles) 
{
	const auto r{ obs_radius }; //Observation point radius distance. We consider a fairly distant enough radius.
	const auto& phi{ angles.phi };
	const auto& theta{ angles.theta };
	const auto x{ r * std::sin(theta) * std::cos(phi) };
	const auto y{ r * std::sin(theta) * std::sin(phi) };
	const auto z{ r * std::cos(theta) };
	const Vector res({ x, y, z });
	return res;
}

double calcPsiAngle3D(const Vector& p, const SphericalAngles& angles)
{
	auto vec{ buildObsPointVec(angles) };
	return std::acos((vec[0] * p[0] + vec[1] * p[1] + vec[2] * p[2]) / (vec.Norml2() * p.Norml2()));
}

double func_exp_real_part_3D(const Vector& p, const Frequency freq, const SphericalAngles& angles)
{
	auto landa = physicalConstants::speedOfLight / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto psi = calcPsiAngle3D(p, angles);
	auto rad_term = wavenumber * p.Norml2() * std::cos(psi);
	return cos(rad_term);
}

double func_exp_imag_part_3D(const Vector& p, const Frequency freq, const SphericalAngles& angles)
{
	auto landa = physicalConstants::speedOfLight / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto psi = calcPsiAngle3D(p, angles);
	auto rad_term = wavenumber * p.Norml2() * std::cos(psi);
	return sin(rad_term);
}

std::unique_ptr<FunctionCoefficient> buildFC_3D(const Frequency freq, const SphericalAngles& angles, bool isReal)
{
	std::function<double(const Vector&)> f;

	if (isReal) {
		f = [freq, angles](const Vector& p) { return func_exp_real_part_3D(p, freq, angles); };
	}
	else {
		f = [freq, angles](const Vector& p) { return func_exp_imag_part_3D(p, freq, angles); };
	}

	FunctionCoefficient res(f);
	return std::make_unique<FunctionCoefficient>(res);
}

Array<int> getNearToFarFieldMarker(const int att_size)
{
	Array<int> res(att_size);
	res = 0;
	res[static_cast<int>(BdrCond::NearToFarField) - 1] = 1;
	return res;
}

std::unique_ptr<LinearForm> assembleLinearForm(ScalarVectorProductCoefficient& vc, FiniteElementSpace& fes)
{
	auto res{ std::make_unique<LinearForm>(&fes) };
	auto marker{ getNearToFarFieldMarker(fes.GetMesh()->bdr_attributes.Max()) };
	res->AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(vc), marker);
	res->Assemble();
	return res;
}

std::unique_ptr<FiniteElementSpace> buildFESFromGF(Mesh& mesh, const std::string& path)
{
	for (auto const& dir_entry : std::filesystem::directory_iterator(path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			std::ifstream inex(dir_entry.path().generic_string() + "/Ex.gf");
			FiniteElementSpace fes;
			fes.Load(&mesh, inex);
			return std::make_unique<FiniteElementSpace>(fes);
		}
	}
	throw std::runtime_error("FES from GridFunciton not returned. Verify input path.");
}

std::map<SphericalAngles, Freq2Value> initAngles2FreqValues(const std::vector<Frequency>& frequencies, const std::vector<SphericalAngles>& angleVec)
{
	std::map<SphericalAngles, Freq2Value> res;
	for (const auto& angpair : angleVec) {
		Freq2Value f2v;
		for (const auto& f : frequencies) {
			f2v.emplace(f, 0.0);
		}
		res.emplace(angpair, f2v);
	}
	return res;
}

GridFunction getGridFunction(Mesh& mesh, const std::string& path)
{
	std::ifstream in(path);
	GridFunction res(&mesh, in);
	return res;
}

const Time getTime(const std::string& timePath)
{
	std::ifstream timeFile(timePath);
	if (!timeFile) {
		throw std::runtime_error("File could not be opened in getTime.");
	}
	std::string timeString;
	std::getline(timeFile, timeString);
	return std::stod(timeString);
}

PlaneWaveData buildPlaneWaveData(const json& json)
{
	double spread(-1e5);
	mfem::Vector mean;
	double projMean(-1e5);

	for (auto s{ 0 }; s < json["sources"].size(); s++) {
		if (json["sources"][s]["type"] == "planewave") {
			spread = json["sources"][s]["magnitude"]["spread"];
			mean = driver::assemble3DVector(json["sources"][s]["magnitude"]["mean"]);
		    mfem::Vector propagation = driver::assemble3DVector(json["sources"][s]["propagation"]);
			projMean = mean * propagation / propagation.Norml2();
		}
	}

	if (std::abs(spread - 1e5) < 1e-6 || std::abs(projMean - 1e5) < 1e-6) {
		throw std::runtime_error("Verify PlaneWaveData inputs for RCS normalization term.");
	}

	return PlaneWaveData(spread, projMean);
}

std::vector<double> buildTimeVector(const std::string& data_path)
{
	std::vector<double> res;
	for (auto const& dir_entry : std::filesystem::directory_iterator(data_path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh" && 
			dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 3) != "rcs"  &&
			dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 8) != "farfield") 
		{
			res.push_back(getTime(dir_entry.path().generic_string() + "/time.txt") / physicalConstants::speedOfLight);
		}
	}
	return res;
}

std::vector<double> evaluateGaussianVector(std::vector<Time>& time, double spread, double mean)
{
	std::vector<double> res(time.size());
	for (int t = 0; t < time.size(); ++t) {
		res[t] = exp(-0.5 * std::pow((time[t] - std::abs(mean)) / spread, 2.0));
	}
	return res;
}

Freq2CompVec calculateDFT(const Vector& gf, const std::vector<Frequency>& frequencies, const Time time)
{
	Freq2CompVec res(frequencies.size());
	for (int f{ 0 }; f < frequencies.size(); f++) {
		res[f].resize(gf.Size(), std::complex<double>(0.0, 0.0));
		for (int i{ 0 }; i < gf.Size(); i++) {
			auto arg = 2.0 * M_PI * frequencies[f] * time;
			auto w = std::complex<double>(cos(arg), -sin(arg));
			res[f][i] += gf[i] * w;
		}
	}
	return res;
}

NedelecFields FarField::buildNedelecFields(Mesh& mesh, FiniteElementSpace& dgfes, const std::string& path)
{

	std::vector<double> time{ buildTimeVector(path) };

	auto dgfec = dynamic_cast<const DG_FECollection*>(dgfes.FEColl());
	auto dgfesv3 = FiniteElementSpace(&mesh, dgfec, 3);
	NedelecFields res(time.size());

	for (auto const& dir_entry : std::filesystem::directory_iterator(path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh" &&
			dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 3) != "rcs" &&
			dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 8) != "farfield") {

				std::ifstream inex(dir_entry.path().generic_string() + "/Ex.gf");
				std::ifstream iney(dir_entry.path().generic_string() + "/Ey.gf");
				std::ifstream inez(dir_entry.path().generic_string() + "/Ez.gf");
				std::ifstream inhx(dir_entry.path().generic_string() + "/Hx.gf");
				std::ifstream inhy(dir_entry.path().generic_string() + "/Hy.gf");
				std::ifstream inhz(dir_entry.path().generic_string() + "/Hz.gf");

				GridFunction gfex(&mesh, inex);
				GridFunction gfey(&mesh, iney);
				GridFunction gfez(&mesh, inez);
				GridFunction gfhx(&mesh, inhx);
				GridFunction gfhy(&mesh, inhy);
				GridFunction gfhz(&mesh, inhz);

				GridFunction dg_gf_electric(&dgfesv3), dg_gf_magnetic(&dgfesv3);

				dg_gf_electric.SetVector(gfex, 0);
				dg_gf_electric.SetVector(gfey, gfex.Size());
				dg_gf_electric.SetVector(gfez, 2 * gfex.Size());
				VectorGridFunctionCoefficient dg_vgfc_electric(&dg_gf_electric);
				GridFunction nd_gf_electric(ndfes_.get());
				nd_gf_electric.ProjectCoefficient(dg_vgfc_electric);
				res.electric.emplace_back(nd_gf_electric);

				dg_gf_magnetic.SetVector(gfhx, 0);
				dg_gf_magnetic.SetVector(gfhy, gfhx.Size());
				dg_gf_magnetic.SetVector(gfhz, 2 * gfhx.Size());
				VectorGridFunctionCoefficient dg_vgfc_magnetic(&dg_gf_magnetic);
				GridFunction nd_gf_magnetic(ndfes_.get());
				nd_gf_magnetic.ProjectCoefficient(dg_vgfc_electric);
				res.magnetic.emplace_back(nd_gf_magnetic);

		}
	}

	return res;
}

RealImagFreqFields FarField::buildFrequencyFields(const NedelecFields& time_fields, const std::vector<Time>& time, const Frequency frequency)
{

	FreqNedelecComponents freqReal, freqImag;

	freqReal.electric.SetSpace(ndfes_.get());
	freqReal.magnetic.SetSpace(ndfes_.get());
	freqImag.electric.SetSpace(ndfes_.get());
	freqImag.magnetic.SetSpace(ndfes_.get());

	freqReal.electric = 0.0;
	freqReal.magnetic = 0.0;
	freqImag.electric = 0.0;
	freqImag.magnetic = 0.0;

	for (int t{ 0 }; t < time.size(); ++t) {
		auto arg = 2.0 * M_PI * frequency * time[t];
		for (int v{ 0 }; v < ndfes_->GetNDofs(); ++v) {
			freqReal.electric[v] += time_fields.electric[t][v] *  cos(arg);
			freqReal.magnetic[v] += time_fields.magnetic[t][v] *  cos(arg);
			freqImag.electric[v] += time_fields.electric[t][v] * -sin(arg);
			freqImag.magnetic[v] += time_fields.magnetic[t][v] * -sin(arg);
		}
	}

	for (int v{ 0 }; v < ndfes_->GetNDofs(); ++v) {
		freqReal.electric[v] /= time.size();
		freqReal.magnetic[v] /= time.size();
		freqImag.electric[v] /= time.size();
		freqImag.magnetic[v] /= time.size();
	}

	return std::pair(freqReal, freqImag);
}

CurrentTerms FarField::integrateCurrents(const RealImagFreqFields& nf, const Frequency frequency, const SphericalAngles& angles)
{
	const auto& nf_real_electric = nf.first.electric;
	const auto& nf_imag_electric = nf.second.electric;
	const auto& nf_real_magnetic = nf.first.magnetic;
	const auto& nf_imag_magnetic = nf.second.magnetic;

	std::unique_ptr<FunctionCoefficient> fc_real, fc_imag;
	switch (dgfes_->GetMesh()->SpaceDimension()) {
	case 3:
		fc_real = buildFC_3D(frequency, angles, true);
		fc_imag = buildFC_3D(frequency, angles, false);
		break;
	default:
		throw std::runtime_error("RCS currently only supported for 3D problems.");
	}

	VectorGridFunctionCoefficient vgfc_e_real(&nf.first.electric), vgfc_e_imag(&nf.second.electric),
								  vgfc_h_real(&nf.first.magnetic), vgfc_h_imag(&nf.second.magnetic);

	GridFunction ones(ndfes_.get()), negones(ndfes_.get());
	ones = 1.0; negones = -1.0;
	VectorGridFunctionCoefficient vcc_ones(&ones), vcc_negones(&negones);

	ScalarVectorProductCoefficient svpc_e_real(*fc_real, vcc_negones);
	ScalarVectorProductCoefficient svpc_e_imag(*fc_imag, vcc_negones);

	ScalarVectorProductCoefficient svpc_h_real(*fc_real, vcc_ones);
	ScalarVectorProductCoefficient svpc_h_imag(*fc_imag, vcc_ones);

	auto lf_M_real = assembleLinearForm(svpc_e_real, *ndfes_);
	auto lf_M_imag = assembleLinearForm(svpc_e_imag, *ndfes_);

	auto lf_J_real = assembleLinearForm(svpc_h_real, *ndfes_);
	auto lf_J_imag = assembleLinearForm(svpc_h_imag, *ndfes_);

	std::complex<double> J_term(lf_J_real.get()->operator()(nf_real_electric) - lf_J_imag.get()->operator()(nf_imag_electric), lf_J_real.get()->operator()(nf_imag_electric) + lf_J_imag.get()->operator()(nf_real_electric));
	std::complex<double> M_term(lf_M_real.get()->operator()(nf_real_magnetic) - lf_M_imag.get()->operator()(nf_imag_magnetic), lf_M_real.get()->operator()(nf_imag_magnetic) + lf_M_imag.get()->operator()(nf_real_magnetic));

	CurrentTerms res;
	res.J = J_term;
	res.M = M_term;
	return res;
}

FarField::FarField(const std::string& data_path, const std::string& json_path, std::vector<Frequency>& frequencies, const std::vector<SphericalAngles>& angle_vec)
{

	auto mesh = Mesh::LoadFromFile(data_path + "/mesh");
	dgfes_ = buildFESFromGF(mesh, data_path);

	auto dgfec = dynamic_cast<const DG_FECollection*>(dgfes_->FEColl());
	auto dgfesv3 = FiniteElementSpace(&mesh, dgfec, 3);
	auto ndfec = ND_FECollection(dgfec->GetOrder(), mesh.Dimension());
	ndfes_ = std::make_unique<FiniteElementSpace>(&mesh, &ndfec);

	pot_rad_ = initAngles2FreqValues(frequencies, angle_vec);

	NedelecFields nedFields{ buildNedelecFields(mesh, *dgfes_ , data_path) };

	auto time = buildTimeVector(data_path);

	Freq2NedFields frequency_fields;
	for (int f = 0; f < frequencies.size(); f++) {
		auto freq_field = buildFrequencyFields(nedFields, time, frequencies[f]);
		frequency_fields[frequencies[f]].first.SetSpace(ndfes_.get());
		frequency_fields[frequencies[f]].second.SetSpace(ndfes_.get());
		frequency_fields[frequencies[f]].first.electric = freq_field.first.electric;
		frequency_fields[frequencies[f]].first.magnetic = freq_field.first.magnetic;
		frequency_fields[frequencies[f]].second.electric = freq_field.second.electric;
		frequency_fields[frequencies[f]].second.magnetic = freq_field.second.magnetic;
	}

	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			auto freq_currents = integrateCurrents(frequency_fields[frequencies[f]], frequencies[f], angpair);
			auto landa = physicalConstants::speedOfLight / frequencies[f];
			auto wavenumber = 2.0 * M_PI / landa;
			auto const_term = std::pow(wavenumber, 2.0) / (4.0 * M_PI);
			auto freq_val = const_term * std::pow(std::abs(physicalConstants::freeSpaceImpedance * freq_currents.J - freq_currents.M),2.0);
			pot_rad_[angpair][frequencies[f]] = freq_val;
		}
	}

}

}
