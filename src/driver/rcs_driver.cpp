#include "rcs_driver.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>

#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "components/FarField.h"
#include "components/RCSSurfacePostProcessor.h"
#include "driver/driver.h"

using json = nlohmann::json;

namespace maxwell::driver {

static std::vector<double> linspace(double a, double b, size_t N)
{
	if (N == 0) return {};
	if (N == 1) return { a };
	std::vector<double> xs(N);
	double h = (b - a) / static_cast<double>(N - 1);
	for (size_t i = 0; i < N; ++i)
		xs[i] = a + i * h;
	return xs;
}

void runRCSPostProcessing(const std::string& rcsJsonPath)
{
	std::ifstream f(rcsJsonPath);
	if (!f)
		throw std::runtime_error("Cannot open RCS input file: " + rcsJsonPath);
	auto rcsInput = json::parse(f);

	std::string runmode  = rcsInput.at("runmode");
	std::string casename = rcsInput.at("casename");

	std::string caseJson = "./testData/maxwellInputs/" + casename + "/" + casename + ".json";

	auto caseData = parseJSONfile(caseJson);

	if (!caseData.contains("probes") || !caseData["probes"].contains("rcssurface"))
		throw std::runtime_error("Case JSON has no rcssurface probes: " + caseJson);

	auto& freqSpec = rcsInput.at("frequencies");
	auto frequencies = linspace(
		freqSpec.at("start").get<double>(),
		freqSpec.at("end").get<double>(),
		freqSpec.at("steps").get<size_t>());

	auto& angSpec = rcsInput.at("angles");
	auto thetas = linspace(
		angSpec.at("theta").at("start").get<double>(),
		angSpec.at("theta").at("end").get<double>(),
		angSpec.at("theta").at("steps").get<size_t>());
	auto phis = linspace(
		angSpec.at("phi").at("start").get<double>(),
		angSpec.at("phi").at("end").get<double>(),
		angSpec.at("phi").at("steps").get<size_t>());

	std::vector<SphericalAngles> angles;
	for (const auto& t : thetas)
		for (const auto& p : phis)
			angles.push_back({ t, p });

	std::optional<double> maxTime;
	if (rcsInput.contains("max_time") && !rcsInput["max_time"].is_null())
		maxTime = rcsInput["max_time"].get<double>();

	for (const auto& probe : caseData["probes"]["rcssurface"]) {
		std::string probeName = probe.at("name");
		std::string dataPath = "./Exports/" + runmode + "/" + casename
			+ "/RCSSurface/" + probeName + "/";

		if (mfem::Mpi::WorldRank() == 0) {
			std::cout << "Processing RCS probe: " << probeName << "\n"
				<< "  Data path: " << dataPath << "\n"
				<< "  Case JSON: " << caseJson << "\n";
		}

		RCSSurfacePostProcessor pp(dataPath, caseJson, frequencies, angles, maxTime);
	}
}

} // namespace maxwell::driver
