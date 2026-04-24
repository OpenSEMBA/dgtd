#include "driver.h"
#include "string"

#include <numeric>
#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <unordered_set>
#include <mpi.h>

namespace maxwell::driver {

double calculateMaximumSourceFrequency(const json& case_data)
{
    double c_si = physicalConstants::speedOfLight_SI;
    double min_spread = std::numeric_limits<double>::max();
    double max_freq = 0.0;
    bool found_gaussian = false;
    bool found_modulated = false;

    if (case_data.contains("sources")) {
        for (const auto& source : case_data["sources"]) {
            if (!source.contains("magnitude")) continue;
            const auto& mag = source["magnitude"];
            // Modulated Gaussian: has "frequency" (with or without explicit "type")
            if (mag.contains("frequency") && mag.contains("spread")) {
                double spread = mag["spread"].get<double>();
                double f_carrier = mag["frequency"].get<double>();
                double f_edge = f_carrier + c_si / (2.0 * spread);
                if (f_edge > max_freq) {
                    max_freq = f_edge;
                    found_modulated = true;
                }
            // Plain Gaussian: explicit type="gaussian", no frequency
            } else if (mag.contains("type") && mag["type"] == "gaussian" && mag.contains("spread")) {
                double spread = mag["spread"].get<double>();
                if (spread > 0.0 && spread < min_spread) {
                    min_spread = spread;
                    found_gaussian = true;
                }
            }
        }
    }

    if (found_modulated) {
        return max_freq;
    }

    if (found_gaussian) {
        return c_si / (2.0 * min_spread);
    }

    return 1e9;
}

std::vector<std::pair<int,int>> buildTwoElementPairsByTagToSort(mfem::Mesh& mesh, mfem::Array<int> tags)
{
	std::vector<std::pair<int,int>> res;
	for (auto t = 0; t < tags.Size(); t++){
		for(auto b = 0; b < mesh.GetNBE(); b++){
			if (mesh.GetBdrAttribute(b) == tags[t]){
				auto f_trans = mesh.GetInternalBdrFaceTransformations(b);
				if (auto f_trans = mesh.GetInternalBdrFaceTransformations(b)) {
    				res.emplace_back(f_trans->Elem1No, f_trans->Elem2No);
				}
			}
		}
	}
	return res;
}

static std::vector<std::pair<int,int>>
gatherConstraintPairs(mfem::Mesh& mesh,
                      const mfem::Array<int>& tfsf_tags,
                      const mfem::Array<int>& sgbc_tags)
{
    std::vector<std::pair<int,int>> pairs;
    pairs.reserve(64);

    if (tfsf_tags.Size() > 0) {
        auto tfsf_pairs = buildTwoElementPairsByTagToSort(mesh, tfsf_tags);
        pairs.insert(pairs.end(), tfsf_pairs.begin(), tfsf_pairs.end());
    }
    if (sgbc_tags.Size() > 0) {
        auto sgbc_pairs = buildTwoElementPairsByTagToSort(mesh, sgbc_tags);
        pairs.insert(pairs.end(), sgbc_pairs.begin(), sgbc_pairs.end());
    }
    return pairs;
}

// Weight for boundary-adjacent elements: SGBC sub-solvers and TFSF source
// injection are significantly more expensive than a plain DG element update.
static constexpr long long BOUNDARY_ELEM_WEIGHT = 5;
static constexpr long long BULK_ELEM_WEIGHT     = 1;

// Build per-element weight array.  Elements adjacent to any SGBC or TFSF
// face get BOUNDARY_ELEM_WEIGHT; all others get BULK_ELEM_WEIGHT.
static std::vector<long long> buildElementWeights(
    int NE,
    const std::vector<std::pair<int,int>>& pairs)
{
    std::vector<long long> w(NE, BULK_ELEM_WEIGHT);
    for (const auto& pr : pairs) {
        if (0 <= pr.first  && pr.first  < NE) w[pr.first]  = BOUNDARY_ELEM_WEIGHT;
        if (0 <= pr.second && pr.second < NE) w[pr.second] = BOUNDARY_ELEM_WEIGHT;
    }
    return w;
}

// Weighted load per rank.
static std::vector<long long> computeWeightedLoad(
    int P, int NE,
    const int* partitioning,
    const std::vector<long long>& weight)
{
    std::vector<long long> load(P, 0);
    for (int e = 0; e < NE; ++e) {
        const int r = partitioning[e];
        if (0 <= r && r < P) load[r] += weight[e];
    }
    return load;
}

// Count how many boundary-pair DOF points (SGBC or TFSF face endpoints)
// are currently assigned to each rank.
static std::vector<int> countBoundaryPairsPerRank(
    int P,
    const int* partitioning,
    const std::vector<std::pair<int,int>>& pairs)
{
    std::vector<int> cnt(P, 0);
    for (const auto& pr : pairs) {
        const int r = partitioning[pr.first];
        if (0 <= r && r < P) cnt[r]++;
    }
    return cnt;
}

// Pick the rank with the lowest weighted load.
static int leastLoadedRank(const std::vector<long long>& load)
{
    int best = 0;
    for (int r = 1; r < static_cast<int>(load.size()); ++r) {
        if (load[r] < load[best]) best = r;
    }
    return best;
}

// Assign a vector of elements to a rank, updating weights.
static void assignComponent(const std::vector<int>& elems,
                             int rank,
                             int* partitioning,
                             std::vector<long long>& load,
                             const std::vector<long long>& weight)
{
    for (int e : elems) {
        const int old_r = partitioning[e];
        if (0 <= old_r && old_r < static_cast<int>(load.size()))
            load[old_r] -= weight[e];
        partitioning[e] = rank;
        load[rank] += weight[e];
    }
}

void applyPairwiseConstraintsPartitioning(mfem::Mesh& mesh,
                                          int* partitioning,
                                          const mfem::Array<int>& tfsf_tags,
                                          const mfem::Array<int>& sgbc_tags)
{
    const int NE = mesh.GetNE();
    const int P  = Mpi::WorldSize();
    if (NE == 0 || P <= 1) return;

    auto pairs = gatherConstraintPairs(mesh, tfsf_tags, sgbc_tags);
    if (pairs.empty()) return;

    // --- Phase 0: Build weighted load model ---
    auto weight = buildElementWeights(NE, pairs);
    auto load   = computeWeightedLoad(P, NE, partitioning, weight);

    // --- Phase 1: Transitive grouping via DSU (union-find) ---
    // At corners a single TF element may share faces with multiple SF
    // elements (or vice-versa).  All such elements must land on the same
    // rank so that every TFSF face is rank-local.  We use a DSU to merge
    // elements transitively through shared boundary faces.
    std::vector<int> parent(NE);
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<int> rank_dsu(NE, 0);

    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (rank_dsu[a] < rank_dsu[b]) std::swap(a, b);
        parent[b] = a;
        if (rank_dsu[a] == rank_dsu[b]) rank_dsu[a]++;
    };

    for (const auto& pr : pairs) {
        const int a = pr.first, b = pr.second;
        if (a < 0 || a >= NE || b < 0 || b >= NE) continue;
        unite(a, b);
    }

    // Collect connected components.
    std::unordered_map<int, std::vector<int>> components;
    for (const auto& pr : pairs) {
        for (int e : {pr.first, pr.second}) {
            if (e >= 0 && e < NE) {
                components[find(e)];  // ensure key exists
            }
        }
    }
    for (auto& [root, elems] : components) {
        elems.clear();
    }
    for (const auto& pr : pairs) {
        for (int e : {pr.first, pr.second}) {
            if (e >= 0 && e < NE) {
                auto& v = components[find(e)];
                if (v.empty() || v.back() != e) {
                    v.push_back(e);
                }
            }
        }
    }
    // Deduplicate element lists (an element may appear in multiple pairs).
    for (auto& [root, elems] : components) {
        std::sort(elems.begin(), elems.end());
        elems.erase(std::unique(elems.begin(), elems.end()), elems.end());
    }

    // Build sortable list of components by descending cost.
    struct Component {
        std::vector<int> elems;
        long long cost;
    };
    std::vector<Component> comps;
    comps.reserve(components.size());
    for (auto& [root, elems] : components) {
        long long c = 0;
        for (int e : elems) c += weight[e];
        comps.push_back({std::move(elems), c});
    }
    std::sort(comps.begin(), comps.end(),
              [](const Component& a, const Component& b) {
                  return a.cost > b.cost;
              });

    // --- Phase 2: Assign each component atomically to least-loaded rank ---
    for (const auto& comp : comps) {
        assignComponent(comp.elems, leastLoadedRank(load),
                        partitioning, load, weight);
    }

    // --- Phase 3: Diagnostics ---
    if (Mpi::WorldRank() == 0) {
        auto bdr_cnt = countBoundaryPairsPerRank(P, partitioning, pairs);
        long long max_load = *std::max_element(load.begin(), load.end());
        long long min_load = *std::min_element(load.begin(), load.end());
        long long total    = std::accumulate(load.begin(), load.end(), 0LL);
        double ideal       = static_cast<double>(total) / P;
        double imbalance   = (ideal > 0.0) ? (max_load - ideal) / ideal * 100.0 : 0.0;

        std::cout << "[Partition] " << pairs.size() << " boundary pairs, "
                  << comps.size() << " components across "
                  << P << " ranks (TFSF=" << tfsf_tags.Size()
                  << " tags, SGBC=" << sgbc_tags.Size() << " tags)\n";
        std::cout << "[Partition] Weighted load: min=" << min_load
                  << " max=" << max_load << " ideal=" << std::fixed
                  << std::setprecision(1) << ideal
                  << " imbalance=" << imbalance << "%\n";
        std::cout << "[Partition] Boundary pairs per rank:";
        for (int r = 0; r < P; ++r) std::cout << " R" << r << "=" << bdr_cnt[r];
        std::cout << std::endl;
    }
}

inline void checkIfThrows(bool condition, const std::string& msg)
{
	if (!condition) {
		throw std::runtime_error(msg.c_str());
	}
}

const FieldType assignFieldType(const std::string& field_type)
{
	if (field_type == "electric") {
		return FieldType::E;
	}
	else if (field_type == "magnetic") {
		return FieldType::H;
	}
	else {
		throw std::runtime_error("Wrong Field Type in Point Probe assignation.");
	}
}

const Direction assignFieldPol(const std::string& direction)
{
	if (direction == "X") {
		return X;
	}
	else if (direction == "Y") {
		return Y;
	}
	else if (direction == "Z") {
		return Z;
	}
	else {
		throw std::runtime_error("Wrong Field Polarization in Point Probe assignation.");
	}
}

std::vector<double> assembleVector(const json& input)
{
	std::vector<double> res(input.size());
	for (int i = 0; i < input.size(); i++) {
		res[i] = input[i];
	}
	return res;
}

mfem::Vector assembleCenterVector(const json& source_center)
{
	mfem::Vector res(int(source_center.size()));
	for (int i = 0; i < source_center.size(); i++) {
		res[i] = source_center[i];
	}
	return res;
}

mfem::Vector assemble3DVector(const json& input)
{
	mfem::Vector res(3);
	for (int i = 0; i < input.size(); i++) {
		res[i] = input[i];
	}
	return res;
}

FieldType getFieldType(const std::string& ft)
{
	if (ft == "electric") {
		return FieldType::E;
	}
	else if (ft == "magnetic") {
		return FieldType::H;
	}
	else {
		throw std::runtime_error("The fieldtype written in the json is neither 'electric' nor 'magnetic'");
	}
}

std::unique_ptr<InitialField> buildGaussianInitialField(
	const FieldType& ft = E,
	const double spread = 0.1,
	const mfem::Vector& center_ = mfem::Vector({ 0.5 }),
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const int dimension = 1)
{
	mfem::Vector gaussianCenter(dimension);
	gaussianCenter = 0.0;

	Gaussian gauss(spread, gaussianCenter, dimension);
	return std::make_unique<InitialField>(gauss, ft, p, center_);
}

std::unique_ptr<InitialField> buildResonantModeInitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const std::vector<std::size_t>& modes = { 1 })
{
	Sources res;
	Source::Position center((int)modes.size());
	center = 0.0;
	return std::make_unique<InitialField>(SinusoidalMode{ modes }, ft, p, center);
}

std::unique_ptr<InitialField> buildBesselJ6InitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0, 0.0, 1.0 }))
{
	Sources res;
	Source::Position center = Source::Position({ 0.0, 0.0, 0.0 });
	return std::make_unique<InitialField>(BesselJ6(), ft, p, center);
}

// ---------------------------------------------------------------------------
// Helpers for automatic source delay (auto-mean) computation
// ---------------------------------------------------------------------------

// How many Gaussian sigma to delay the source past the TFSF arrival time.
// exp(-N^2) ~ 1.4e-11 for N=5, ensuring negligible field at t=0.
static constexpr double AUTO_DELAY_N_SIGMA = 5.0;

// Returns the centroid of boundary element `be` using its vertex coordinates.
static mfem::Vector bdrElemCentroid(const mfem::Mesh& mesh, int be)
{
	mfem::Array<int> verts;
	mesh.GetBdrElementVertices(be, verts);
	int dim = mesh.Dimension();
	mfem::Vector c(dim);
	c = 0.0;
	for (int v = 0; v < verts.Size(); ++v) {
		const double* p = mesh.GetVertex(verts[v]);
		for (int d = 0; d < dim; ++d) c[d] += p[d];
	}
	c /= static_cast<double>(verts.Size());
	return c;
}

// Returns true if boundary element `be` has one of the given JSON tags.
static bool bdrElemHasTag(const mfem::Mesh& mesh, int be, const json& tags)
{
	int attr = mesh.GetBdrAttribute(be);
	for (const auto& t : tags) {
		if (t.get<int>() == attr) return true;
	}
	return false;
}

// Minimum Euclidean distance from the origin of all TFSF surface elements.
// For dipole sources: the retarded-time wave must travel at least this far
// before it reaches the TFSF surface, so it constrains the required mean.
static double minRadiusOnTFSFSurface(const mfem::Mesh& mesh, const json& tags)
{
	double min_r = std::numeric_limits<double>::max();
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		if (!bdrElemHasTag(mesh, be, tags)) continue;
		min_r = std::min(min_r, bdrElemCentroid(mesh, be).Norml2());
	}
	double global_min_r;
	MPI_Allreduce(&min_r, &global_min_r, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return (global_min_r == std::numeric_limits<double>::max()) ? 0.0 : global_min_r;
}

// Minimum phase delay x·d̂/c of all TFSF surface elements in propagation direction d_hat.
// For a planewave traveling in direction d_hat, this is the phase at the
// earliest-arriving (most upstream) point of the TFSF surface.  The auto-mean
// is chosen so that the pulse is N sigma before peak at this earliest point.
static double minPhaseOnTFSFSurface(const mfem::Mesh& mesh, const json& tags,
	const mfem::Vector& d_hat)
{
	double min_phase = std::numeric_limits<double>::max();
	const double c = physicalConstants::speedOfLight;
	int dim = mesh.Dimension();
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		if (!bdrElemHasTag(mesh, be, tags)) continue;
		auto x = bdrElemCentroid(mesh, be);
		double phase = 0.0;
		for (int d = 0; d < dim && d < d_hat.Size(); ++d) phase += x[d] * d_hat[d];
		phase /= c;
		min_phase = std::min(min_phase, phase);
	}
	double global_min_phase;
	MPI_Allreduce(&min_phase, &global_min_phase, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return (global_min_phase == std::numeric_limits<double>::max()) ? 0.0 : global_min_phase;
}

std::unique_ptr<InitialField> buildSphericalBesselJ6InitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0, 0.0, 1.0 }))
{
	Sources res;
	Source::Position center = Source::Position({ 0.0, 0.0, 0.0 });
	return std::make_unique<InitialField>(SphericalBesselJ6(), ft, p, center);
}

std::unique_ptr<TotalField> buildGaussianPlanewave(
	double spread,
	const Source::Position mean,
	const Source::Polarization& pol,
	const Source::Propagation& dir,
	const FieldType ft = FieldType::E
)
{
	Position projMean(3);
	projMean = 0.0;
	for (auto v = 0; v < mean.Size(); v++) {
		projMean[v] = mean[v];
	}
	Gaussian gauss{ spread, mfem::Vector({projMean * dir / dir.Norml2()})};
	Planewave pw(gauss, pol, dir, ft);
	return std::make_unique<TotalField>(pw);
}

std::unique_ptr<TotalField> buildModulatedGaussianPlanewave(
	double spread,
	const Source::Position mean,
	double freq_hz,
	const Source::Polarization& pol,
	const Source::Propagation& dir,
	const FieldType ft = FieldType::E
)
{
	Position projMean(3);
	projMean = 0.0;
	for (auto v = 0; v < mean.Size(); v++) {
		projMean[v] = mean[v];
	}
	double freq_norm = freq_hz / physicalConstants::speedOfLight_SI;
	ModulatedGaussian mg{ spread, mfem::Vector({projMean * dir / dir.Norml2()}), freq_norm };
	Planewave pw(mg, pol, dir, ft);
	return std::make_unique<TotalField>(pw);
}

std::unique_ptr<TotalField> buildDerivGaussDipole(
	const double length, 
	const double gaussianSpread, 
	const double gaussMean
) 
{
	DerivGaussDipole dip(length, gaussianSpread, gaussMean);
	return std::make_unique<TotalField>(dip);
}

Sources buildSources(const json& case_data, const mfem::Mesh* mesh)
{
	Sources res;
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		if (case_data["sources"][s]["type"] == "initial") {
			if (case_data["sources"][s]["magnitude"]["type"] == "gaussian") {
				res.add(buildGaussianInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					case_data["sources"][s]["magnitude"]["spread"],
					assembleCenterVector(case_data["sources"][s]["center"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["dimension"])
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "resonant") {
				res.add(buildResonantModeInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["magnitude"]["modes"])
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "besselj6_2D") {
				res.add(buildBesselJ6InitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]))
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "besselj6_3D") {
				res.add(buildSphericalBesselJ6InitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]))
				);
			}
		}
		else if (case_data["sources"][s]["type"] == "planewave") {
			const auto& mag = case_data["sources"][s]["magnitude"];
			double spread = mag["spread"].get<double>();

			// Determine mean: use explicit value if provided, otherwise auto-compute
			// from the earliest-arriving (upstream) point of the TFSF surface so that
			// the pulse is AUTO_DELAY_N_SIGMA sigma below peak at t=0 everywhere.
			Source::Position mean_vec(3);
			mean_vec = 0.0;
			if (mag.contains("mean")) {
				mean_vec = assemble3DVector(mag["mean"]);
			} else if (mesh && case_data["sources"][s].contains("tags")) {
				auto d_hat = assemble3DVector(case_data["sources"][s]["propagation"]);
				d_hat /= d_hat.Norml2();
				double min_phase = minPhaseOnTFSFSurface(*mesh, case_data["sources"][s]["tags"], d_hat);
				// mean_1D = min_phase - N*sigma*sqrt(2): places the pulse N sigma
				// before reaching the upstream TFSF surface at t=0.
				double mean_1D = min_phase - AUTO_DELAY_N_SIGMA * spread * std::sqrt(2.0);
				for (int d = 0; d < d_hat.Size(); ++d) mean_vec[d] = mean_1D * d_hat[d];
				int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				if (rank == 0) {
					std::cout << "[Source " << s << " (planewave)] Auto-computed mean_1D = "
					          << mean_1D << " (min_phase=" << min_phase << ")\n";
				}
			}

			if (mag.contains("frequency")) {
				res.add(buildModulatedGaussianPlanewave(
					spread,
					mean_vec,
					mag["frequency"].get<double>(),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					assemble3DVector(case_data["sources"][s]["propagation"]),
					FieldType::E)
				);
			} else {
				res.add(buildGaussianPlanewave(
					spread,
					mean_vec,
					assemble3DVector(case_data["sources"][s]["polarization"]),
					assemble3DVector(case_data["sources"][s]["propagation"]),
					FieldType::E)
				);
			}
		}
		else if (case_data["sources"][s]["type"] == "dipole") {
			const auto& mag = case_data["sources"][s]["magnitude"];
			double length = mag["length"].get<double>();
			double spread = mag["spread"].get<double>();

			// Determine mean: use explicit value if provided, otherwise auto-compute.
			// For the retarded-time Gaussian, mean_auto ensures the field is
			// AUTO_DELAY_N_SIGMA sigma before peak at t=0 on the TFSF surface.
			double mean;
			if (mag.contains("mean")) {
				mean = mag["mean"].get<double>();
			} else if (mesh && case_data["sources"][s].contains("tags")) {
				double min_r = minRadiusOnTFSFSurface(*mesh, case_data["sources"][s]["tags"]);
				// mean_auto = N*sigma*sqrt(2) - r_min/c
				// (r_min/c is the retarded-time advance the TFSF surface already provides)
				double mean_auto = AUTO_DELAY_N_SIGMA * spread * std::sqrt(2.0)
				                   - min_r / physicalConstants::speedOfLight;
				mean = std::max(spread * std::sqrt(2.0), mean_auto);
				std::cout << "[Source " << s << " (dipole)] Auto-computed mean = "
				          << mean << " (min_radius=" << min_r << ")\n";
			} else {
				mean = AUTO_DELAY_N_SIGMA * spread * std::sqrt(2.0);
			}

			res.add(buildDerivGaussDipole(length, spread, mean));
		}
		else {
			throw std::runtime_error("Unknown source type in Json.");
		}
	}
	return res;
}


SolverOptions buildSolverOptions(const json& case_data)
{
	
	SolverOptions res{};

	if (case_data.contains("solver_options")) {

		if (case_data["solver_options"].contains("upwind_alpha")) {
			res.setUpwindAlpha(case_data["solver_options"]["upwind_alpha"]);
		}

		if (case_data["solver_options"].contains("time_step")) {
			res.setTimeStep(case_data["solver_options"]["time_step"]);
		}

		if (case_data["solver_options"].contains("final_time")) {
			res.setFinalTime(case_data["solver_options"]["final_time"]);
		}

		if (case_data["solver_options"].contains("cfl")) {
			res.setCFL(case_data["solver_options"]["cfl"]);
		}

		if (case_data["solver_options"].contains("order")) {
			res.setOrder(case_data["solver_options"]["order"]);
		}

		if (case_data["solver_options"].contains("spectral")) {
			res.setSpectralEO(case_data["solver_options"]["spectral"]);
		}
		
		if (case_data["solver_options"].contains("export_operator")) {
			res.setExportEO(case_data["solver_options"]["export_operator"]);
		}

		if (case_data["solver_options"].contains("evolution_operator")) {
			if (case_data["solver_options"]["evolution_operator"] == "maxwell") {
				res.setEvolutionOperator(EvolutionOperatorType::Maxwell);
			}
			else if (case_data["solver_options"]["evolution_operator"] == "global") {
				res.setEvolutionOperator(EvolutionOperatorType::Global);
			}
			else if (case_data["solver_options"]["evolution_operator"] == "hesthaven") {
				res.setEvolutionOperator(EvolutionOperatorType::Hesthaven);
			}
			else {
				throw std::runtime_error("Wrong type of Evolution Operator defined, please choose: 'maxwell', 'global' or 'hesthaven'. If none explicitly stated, it will default to 'global'.");
			}
		}

		if (case_data["solver_options"].contains("basis_type")) {
			res.setBasisType(case_data["solver_options"]["basis_type"]);
		}

		if (case_data["solver_options"].contains("ode_type")){
			res.setODEType(case_data["solver_options"]["ode_type"]);
		}

	}

	return res;
}

Probes buildProbes(const json& case_data)
{
    Probes probes;
    size_t pp_count = 0;
    size_t fp_count = 0;

    // Smart step calculator lambda.
    // When time_step == 0.0 the solver will use an automatic time step that is
    // not yet known here, so we return 1 as a safe placeholder; the solver will
    // call ProbesManager::recalculateExportSteps() once the real dt is known.
    auto calculate_interval = [&](const json& probe_data) -> int {
        if (probe_data.contains("steps")) {
            return probe_data["steps"]; // Legacy manual step interval
        }
        else if (probe_data.contains("saves")) {
            int requested_saves = probe_data["saves"];
            if (requested_saves <= 0) return 1;

            double t_final = case_data["solver_options"]["final_time"].get<double>();
            double dt = case_data["solver_options"]["time_step"].get<double>();
            if (dt <= 0.0) return 1; // Placeholder; will be corrected by recalculateExportSteps

            int total_simulation_steps = static_cast<int>(std::ceil(t_final / dt));
            return std::max(1, total_simulation_steps / requested_saves);
        }
        return 1; // Default to every step
    };

    if (case_data.contains("probes")){
        if (case_data["probes"].contains("exporter")) {
            ExporterProbe exporter_probe;
            if (case_data["probes"]["exporter"].contains("name")) {
                exporter_probe.name = case_data["probes"]["exporter"]["name"];
            } else {
                exporter_probe.name = case_data["model"]["filename"];
            }
            exporter_probe.visSteps = calculate_interval(case_data["probes"]["exporter"]);
            if (case_data["probes"]["exporter"].contains("saves"))
                exporter_probe.saves = case_data["probes"]["exporter"]["saves"];
            probes.exporterProbes.push_back(exporter_probe);
        }

        if (case_data["probes"].contains("point")) {
            for (int p = 0; p < case_data["probes"]["point"].size(); p++) {
                int interval = calculate_interval(case_data["probes"]["point"][p]);
                PointProbe point_probe(
                    assembleVector(case_data["probes"]["point"][p]["position"]),
                    interval
                );
                point_probe.setProbeID(pp_count);
                if (case_data["probes"]["point"][p].contains("saves"))
                    point_probe.setSaves(case_data["probes"]["point"][p]["saves"]);
                pp_count++;
                probes.pointProbes.push_back(point_probe);
            }
        }

        if (case_data["probes"].contains("field")) {
            for (int p = 0; p < case_data["probes"]["field"].size(); p++) {
                int interval = calculate_interval(case_data["probes"]["field"][p]);
                FieldProbe field_probe(
                    assignFieldType(case_data["probes"]["field"][p]["field_type"]),
                    assignFieldPol(case_data["probes"]["field"][p]["polarization"]),
                    assembleVector(case_data["probes"]["field"][p]["position"]),
                    interval
                );
                field_probe.setProbeID(fp_count);
                if (case_data["probes"]["field"][p].contains("saves"))
                    field_probe.setSaves(case_data["probes"]["field"][p]["saves"]);
                fp_count++;
                probes.fieldProbes.push_back(field_probe);
            }
        }

        if (case_data["probes"].contains("farfield")) {
            for (int p = 0; p < case_data["probes"]["farfield"].size(); p++) {
                NearFieldProbe probe;
                if (case_data["probes"]["farfield"][p].contains("name")) {
                    probe.name = case_data["probes"]["farfield"][p]["name"];
                }
                if (case_data["probes"]["farfield"][p].contains("export_path")) {
                    probe.exportPath = case_data["probes"]["farfield"][p]["export_path"];
                }
                probe.expSteps = calculate_interval(case_data["probes"]["farfield"][p]);
                if (case_data["probes"]["farfield"][p].contains("saves"))
                    probe.saves = case_data["probes"]["farfield"][p]["saves"];
                
                if (case_data["probes"]["farfield"][p].contains("tags")) {
                    std::vector<int> tags;
                    for (int t = 0; t < case_data["probes"]["farfield"][p]["tags"].size(); t++) {
                        tags.push_back(case_data["probes"]["farfield"][p]["tags"][t]);
                    }
                    probe.tags = tags;
                }
                else {
                    throw std::runtime_error("Tags have not been defined in farfield probe.");
                }
                probes.nearFieldProbes.push_back(probe);
            }
        }

        if (case_data["probes"].contains("domain_snapshot")){
            DomainSnapshotProbe probe;
            if (case_data["probes"]["domain_snapshot"].contains("name")) {
                probe.name = case_data["probes"]["domain_snapshot"]["name"];
            }
            else{
                probe.name = case_data["model"]["filename"];
            }
            probe.expSteps = calculate_interval(case_data["probes"]["domain_snapshot"]);
            if (case_data["probes"]["domain_snapshot"].contains("saves"))
                probe.saves = case_data["probes"]["domain_snapshot"]["saves"];
            probes.domainSnapshotProbes.push_back(probe);
        }

        if (case_data["probes"].contains("rcssurface")) {
            for (size_t p = 0; p < case_data["probes"]["rcssurface"].size(); p++) {
                RCSSurfaceProbe probe;
                if (case_data["probes"]["rcssurface"][p].contains("name")) {
                    probe.name = case_data["probes"]["rcssurface"][p]["name"];
                }
                probe.expSteps = calculate_interval(case_data["probes"]["rcssurface"][p]);
                if (case_data["probes"]["rcssurface"][p].contains("saves"))
                    probe.saves = case_data["probes"]["rcssurface"][p]["saves"];
                if (case_data["probes"]["rcssurface"][p].contains("tags")) {
                    for (size_t t = 0; t < case_data["probes"]["rcssurface"][p]["tags"].size(); t++) {
                        probe.tags.push_back(case_data["probes"]["rcssurface"][p]["tags"][t]);
                    }
                } else {
                    throw std::runtime_error("Tags have not been defined in rcssurface probe.");
                }
                probes.rcsSurfaceProbes.push_back(probe);
            }
        }

        if (case_data["probes"].contains("mor_state")) {
            for (size_t p = 0; p < case_data["probes"]["mor_state"].size(); p++) {
                MORStateProbe probe;
                const auto& pd = case_data["probes"]["mor_state"][p];
                if (pd.contains("name")) {
                    probe.name = pd["name"];
                }
                if (pd.contains("record_time_start")) {
                    probe.record_time_start = pd["record_time_start"].get<double>();
                }
                if (pd.contains("record_time_final")) {
                    probe.record_time_final = pd["record_time_final"].get<double>();
                }
                if (pd.contains("saves")) {
                    probe.saves = pd["saves"].get<int>();
                }
                probes.morStateProbes.push_back(probe);
            }
        }
    }

    return probes;
}


std::string dataFolder() { return "./testData/"; }
std::string maxwellInputsFolder() { return dataFolder() + "maxwellInputs/"; }

BdrCond assignBdrCond(const std::string& bdr_cond)
{
	if (bdr_cond == "PEC") {
		return BdrCond::PEC;
	}
	else if (bdr_cond == "PMC") {
		return BdrCond::PMC;
	}
	else if (bdr_cond == "SMA") {
		return BdrCond::SMA;
	}
	else if (bdr_cond == "SGBC") {
		return BdrCond::SGBC;
	}
	else {
		throw std::runtime_error(("The defined Boundary Type " + bdr_cond + " is incorrect.").c_str());
	}
}

std::string assembleMeshString(const std::string& filename)
{
	std::string folder_name{ filename };
	std::string s_msh = ".msh";
	std::string s_mesh = ".mesh";

	std::string::size_type input_msh = folder_name.find(s_msh);
	std::string::size_type input_mesh = folder_name.find(s_mesh);

	if (input_msh != std::string::npos)
		folder_name.erase(input_msh, s_msh.length());
	if (input_mesh != std::string::npos)
		folder_name.erase(input_mesh, s_mesh.length());

	return maxwellInputsFolder() + folder_name + "/" + filename;
}

std::string assembleLauncherMeshString(const std::string& mesh_name, const std::string& case_path)
{
	std::string path{ case_path };
	std::string s_json = ".json";
	std::string s_msh  = ".msh";
	std::string s_mesh = ".mesh";

	std::string::size_type input_json = path.find(s_json);
	std::string::size_type input_msh = mesh_name.find(s_msh);
	std::string::size_type input_mesh = mesh_name.find(s_mesh);

	path.erase(input_json, s_json.length());
	if (input_msh != std::string::npos)
		path.append(s_msh);
	if (input_mesh != std::string::npos)
		path.append(s_mesh);

	return path;
}

void checkIfAttributesArePresent(const Mesh& mesh, const GeomTagToMaterialInfo& info)
{

	for (auto [att, v] : info.gt2m) {
		checkIfThrows(
			mesh.attributes.Find(att) != -1,
			std::string("There is no attribute") + std::to_string(att) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}

	for (auto [bdr_att, v] : info.gt2bm) {
		checkIfThrows(
			mesh.bdr_attributes.Find(bdr_att) != -1,
			std::string("There is no bdr_attribute") + std::to_string(bdr_att) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}
}

GeomTagToMaterialInfo assembleAttributeToMaterial(const json& case_data, const mfem::Mesh& mesh)
{
	GeomTagToMaterialInfo res{};

	checkIfThrows(case_data.contains("model"), "JSON data does not include 'model'.");
	checkIfThrows(case_data["model"].contains("materials"), "JSON data does not include 'materials'.");

	for (auto m = 0; m < case_data["model"]["materials"].size(); m++) {
		for (auto t = 0; t < case_data["model"]["materials"][m]["tags"].size(); t++) {
			double eps{ 1.0 }, mu{ 1.0 }, sigma{ 0.0 };
			if (case_data["model"]["materials"][m].contains("relative_permittivity")) {
				eps = case_data["model"]["materials"][m]["relative_permittivity"];
			}
			if (case_data["model"]["materials"][m].contains("relative_permeability")) {
				mu = case_data["model"]["materials"][m]["relative_permeability"];
			}
			if (case_data["model"]["materials"][m].contains("bulk_conductivity")) {
				sigma = case_data["model"]["materials"][m]["bulk_conductivity"].get<double>() * physicalConstants::freeSpaceImpedance_SI;
			}
			res.gt2m.emplace(std::make_pair(case_data["model"]["materials"][m]["tags"][t], Material(eps, mu, sigma)));
		}
	}

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			if (case_data["model"]["boundaries"][b].contains("material")) {
				double eps{ 1.0 }, mu{ 1.0 }, sigma{ 0.0 };
				if (case_data["model"]["boundaries"][b]["material"].contains("relative_permittivity")) {
					eps = case_data["model"]["boundaries"][b]["material"]["relative_permittivity"];
				}
				if (case_data["model"]["boundaries"][b]["material"].contains("relative_permeability")) {
					mu = case_data["model"]["boundaries"][b]["material"]["relative_permeability"];
				}
				if (case_data["model"]["boundaries"][b]["material"].contains("bulk_conductivity")) {
					sigma = case_data["model"]["boundaries"][b]["material"]["bulk_conductivity"].get<double>() * physicalConstants::freeSpaceImpedance_SI;
				}
				res.gt2bm.emplace(case_data["model"]["boundaries"][b]["tags"][a], Material(eps, mu, sigma));
			}
		}
	}

	checkIfAttributesArePresent(mesh, res);

	return res;
}

void checkBoundaryInputProperties(const json& case_data)
{
	checkIfThrows(case_data["model"].contains("boundaries"),
		"JSON data does not include 'boundaries' in 'model'.");

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {

		checkIfThrows(case_data["model"]["boundaries"][b].contains("tags"),
			"Boundary " + std::to_string(b) + " does not have defined 'tags'.");

		checkIfThrows(!case_data["model"]["boundaries"][b]["tags"].empty(),
			"Boundary " + std::to_string(b) + " 'tags' are empty.");

		checkIfThrows(case_data["model"]["boundaries"][b].contains("type"),
			"Boundary " + std::to_string(b) + " does not have a defined 'type'.");
	}
}

GeomTagToBoundaryInfo assembleAttributeToBoundary(const json& case_data, const mfem::Mesh& mesh)
{
	using isInterior = bool;

	struct geomTag2Info {
		std::map<GeomTag, BdrCond> geomTag2BdrCond;
		std::map<GeomTag, isInterior> geomTag2IsInterior;
	};

	checkBoundaryInputProperties(case_data);
	auto face2BdrEl{ mesh.GetFaceToBdrElMap() };

	geomTag2Info gt2i;

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			gt2i.geomTag2BdrCond.emplace(case_data["model"]["boundaries"][b]["tags"][a], assignBdrCond(case_data["model"]["boundaries"][b]["type"]));
			gt2i.geomTag2IsInterior.emplace(case_data["model"]["boundaries"][b]["tags"][a], false);
			for (auto f = 0; f < mesh.GetNumFaces(); f++) {
				if (face2BdrEl[f] != -1) {
					if (mesh.GetBdrAttribute(face2BdrEl[f]) == case_data["model"]["boundaries"][b]["tags"][a] 
						&& mesh.FaceIsInterior(f)) {
						gt2i.geomTag2IsInterior[case_data["model"]["boundaries"][b]["tags"][a]] = true;
						break;
					}
				}
			}
		}
	}

	GeomTagToBoundaryInfo res;
	for (auto [att, isInt] : gt2i.geomTag2IsInterior) {
		switch (isInt) {
		case false:
			res.gt2b.emplace(att, gt2i.geomTag2BdrCond[att]);
			break;
		case true:
			res.gt2ib.emplace(att, gt2i.geomTag2BdrCond[att]);
			break;
		}
	}

	return res;
}

// Gmsh 2.2 element type IDs for lower-dimensional entities to strip.
// 1D: nothing stripped.  2D: points (15).  3D: points (15) + lines (1, 8, 26, 27, 28).
static const std::unordered_set<int> kPointTypes   = {15};
static const std::unordered_set<int> kLineTypes     = {1, 8, 26, 27, 28};
// Element types that define the problem dimension.
static const std::unordered_set<int> kTetTypes      = {4, 11, 29};  // tet4, tet10, tet20
static const std::unordered_set<int> kTriTypes      = {2, 9, 20, 21}; // tri3, tri6, tri9, tri10

static int detectMeshDimension(const std::vector<std::string>& element_lines)
{
    int dim = 1;
    for (const auto& line : element_lines) {
        std::istringstream iss(line);
        int id, etype;
        if (!(iss >> id >> etype)) continue;
        if (kTetTypes.count(etype))  return 3;
        if (kTriTypes.count(etype))  dim = std::max(dim, 2);
    }
    return dim;
}

static bool shouldDeleteElement(int etype, int dim)
{
    if (dim == 3) return kPointTypes.count(etype) || kLineTypes.count(etype);
    if (dim == 2) return kPointTypes.count(etype);
    return false; // 1D: keep everything
}

void fixGmshMesh(const std::string& filepath)
{
    // Only process .msh files (Gmsh 2.2 ASCII).
    {
        auto dot = filepath.rfind('.');
        if (dot == std::string::npos || filepath.substr(dot) != ".msh")
            return;
    }

    std::vector<std::string> header_lines;   // everything before first element line
    std::vector<std::string> element_lines;  // raw element lines
    std::vector<std::string> footer_lines;   // $EndElements onward

    {
        std::ifstream in(filepath);
        if (!in.is_open()) return;

        std::string line;
        bool in_elements = false;
        bool reading_elements = false;
        int element_count_line_idx = -1;

        while (std::getline(in, line)) {
            if (line.find("$Elements") != std::string::npos && !in_elements) {
                in_elements = true;
                header_lines.push_back(line);
                // Next line is the element count.
                if (!std::getline(in, line)) break;
                header_lines.push_back(line); // placeholder — will be rewritten
                element_count_line_idx = static_cast<int>(header_lines.size()) - 1;
                reading_elements = true;
                continue;
            }
            if (line.find("$EndElements") != std::string::npos) {
                reading_elements = false;
                footer_lines.push_back(line);
                continue;
            }
            if (reading_elements) {
                element_lines.push_back(line);
            } else if (footer_lines.empty()) {
                header_lines.push_back(line);
            } else {
                footer_lines.push_back(line);
            }
        }
    }

    if (element_lines.empty()) return;

    int dim = detectMeshDimension(element_lines);

    // Process elements: filter + swap physical tag with elementary tag + renumber.
    std::vector<std::string> fixed_elements;
    fixed_elements.reserve(element_lines.size());
    int new_id = 1;

    for (const auto& raw : element_lines) {
        std::vector<std::string> tokens;
        std::istringstream iss(raw);
        std::string tok;
        while (iss >> tok) tokens.push_back(tok);
        if (tokens.size() < 5) continue;

        int etype = std::stoi(tokens[1]);
        if (shouldDeleteElement(etype, dim)) continue;

        // Swap: tokens[3] = physical tag, tokens[4] = elementary tag.
        tokens[3] = tokens[4];
        tokens[0] = std::to_string(new_id++);

        std::ostringstream oss;
        for (size_t i = 0; i < tokens.size(); ++i) {
            if (i) oss << ' ';
            oss << tokens[i];
        }
        fixed_elements.push_back(oss.str());
    }

    // Update element count in header.
    // The element count line is the last line pushed before elements started.
    auto& count_line = header_lines.back();
    count_line = std::to_string(fixed_elements.size());

    // Write back.
    {
        std::ofstream out(filepath, std::ios::trunc);
        for (const auto& l : header_lines) out << l << '\n';
        for (const auto& l : fixed_elements) out << l << '\n';
        for (const auto& l : footer_lines) out << l << '\n';
    }
}

mfem::Mesh assembleMesh(const std::string& mesh_string)
{
	if (Mpi::WorldRank() == 0) {
		fixGmshMesh(mesh_string);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return mfem::Mesh::LoadFromFile(mesh_string, 1, 0, true);
}

mfem::Mesh assembleMeshNoFix(const std::string& mesh_string)
{
	return mfem::Mesh::LoadFromFileNoBdrFix(mesh_string, 1, 0, false);
}

Array<int> getTFSFTags(const json& case_data)
{
    Array<int> res;
    if (!case_data.contains("sources")) return res;

    for (int s = 0; s < case_data["sources"].size(); ++s) {
        if (case_data["sources"][s].contains("type") &&
            (case_data["sources"][s]["type"] == "planewave" ||
             case_data["sources"][s]["type"] == "dipole")) {
            for (int t = 0; t < case_data["sources"][s]["tags"].size(); ++t) {
                res.Append(case_data["sources"][s]["tags"][t]);
            }
        }
    }
    return res;
}

Array<int> getSGBCTags(const json& case_data)
{
    Array<int> res;
    if (!case_data.contains("model") ||
        !case_data["model"].contains("boundaries")) {
        return res;
    }

    for (int b = 0; b < case_data["model"]["boundaries"].size(); ++b) {
        if (case_data["model"]["boundaries"][b].contains("type") && case_data["model"]["boundaries"][b]["type"] == "SGBC") {
            for (int t = 0; t < case_data["model"]["boundaries"][b]["tags"].size(); ++t) {
                res.Append(case_data["model"]["boundaries"][b]["tags"][t]);
            }
        }
    }
    return res;
}

BdrCond assignBoundaryType(const std::string& sgbc_bdr_type)
{
	if (sgbc_bdr_type == "PEC"){
		return BdrCond::PEC;
	}
	if (sgbc_bdr_type == "PMC"){
		return BdrCond::PMC;
	}
	if (sgbc_bdr_type == "SMA"){
		return BdrCond::SMA;
	}
	else{
		throw std::runtime_error("Incorrect sgbc_bdr_type defined in .json.");
	}

}

Model buildModel(const json& case_data, const std::string& case_path, const bool isTest)
{
    mfem::Mesh mesh;
    if (isTest) {
        mesh = assembleMesh(assembleMeshString(case_data["model"]["filename"]));
    } else {
        mesh = assembleMesh(assembleLauncherMeshString(case_data["model"]["filename"], case_path));
    }

    auto att_to_material{ assembleAttributeToMaterial(case_data, mesh) };
    auto att_to_bdr_info{ assembleAttributeToBoundary(case_data, mesh) };

    if (!att_to_bdr_info.gt2b.empty() && !att_to_material.gt2m.empty()) {
        if (isTest) {
            mesh = assembleMesh(assembleMeshString(case_data["model"]["filename"]));
        } else {
            mesh = assembleMesh(assembleLauncherMeshString(case_data["model"]["filename"], case_path));
        }
    }

    if (case_data["model"].contains("refinement")) {
        auto ref_levels = int(case_data["model"]["refinement"]);
        for (auto r = 0; r < ref_levels; r++) {
            mesh.UniformRefinement();
        }
    }

    int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
    mfem::Array<int> tfsf_tags = getTFSFTags(case_data);
    mfem::Array<int> sgbc_tags  = getSGBCTags(case_data);

    if (tfsf_tags.Size() || sgbc_tags.Size()){
        applyPairwiseConstraintsPartitioning(mesh, partitioning, tfsf_tags, sgbc_tags);
    }

    Model res(mesh, att_to_material, att_to_bdr_info, partitioning);
    std::string filename = case_data["model"]["filename"];

    auto ends_with = [](const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() &&
               str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    if (ends_with(filename, ".msh")) {
        res.meshName_ = filename.substr(0, filename.size() - 4);
    } else if (ends_with(filename, ".mesh")) {
        res.meshName_ = filename.substr(0, filename.size() - 5);
    } else {
        throw std::runtime_error("File format for mesh must be name.msh or name.mesh");
    }

    std::vector<SGBCProperties> sgbc_props;
    std::vector<std::string> sgbc_notices;
    
    double max_freq = calculateMaximumSourceFrequency(case_data);

    auto parseSGBCLayer = [&](const nlohmann::json& mat_json) -> SGBCLayer {
        double rel_eps = 1.0;
        if (mat_json.contains("relative_permittivity")) {
            rel_eps = mat_json["relative_permittivity"].get<double>();
        } else {
            std::cout << "SGBC layer defined without 'relative_permittivity', assuming vacuum." << std::endl;
        }

        double rel_mu = 1.0;
        if (mat_json.contains("relative_permeability")) {
            rel_mu = mat_json["relative_permeability"].get<double>();
        } else {
            std::cout << "SGBC layer defined without 'relative_permeability', assuming vacuum." << std::endl;
        }

        if (!mat_json.contains("bulk_conductivity")) {
            throw std::runtime_error("SGBC layer defined without 'bulk_conductivity' parameter. Verify .json parameters.");
        }
        double sigma_si = mat_json["bulk_conductivity"].get<double>();
        double sigma_solver = sigma_si * physicalConstants::freeSpaceImpedance_SI;

        Material mat(rel_eps, rel_mu, sigma_solver);

        if (!mat_json.contains("material_width")) {
            throw std::runtime_error("SGBC layer must define 'material_width'.");
        }
        double layer_width = mat_json["material_width"].get<double>();

        SGBCLayer layer(mat, layer_width);

        double mu_si = rel_mu * physicalConstants::vacuumPermeability_SI;
        double eps_si = rel_eps * physicalConstants::vacuumPermittivity_SI;

        double peclet_number = 0.0;
        if (sigma_si > 0.0) {
            peclet_number = sigma_si / (max_freq * eps_si);
        }

        // Adaptive polynomial order per layer
        if (peclet_number > 100.0) {
            layer.order = 2;
        } else if (peclet_number < 0.1 && sigma_si < 1e-3) {
            layer.order = 4;
        } else {
            layer.order = 3;
        }

        // Segment count per layer
        // Target 20 DOFs per wavelength for accurate phase representation.
        // Each segment has 'order' DOFs, so need ceil(20/order) segments per λ.
        double wavelength = 1.0 / (max_freq * std::sqrt(mu_si * eps_si));
        int segs_per_wavelength = static_cast<int>(std::ceil(20.0 / layer.order));
        double target_dx_wave = wavelength / segs_per_wavelength;

        double target_dx_skin = std::numeric_limits<double>::max();
        double skin_depth = 0.0;
        if (sigma_si > 0.0) {
            skin_depth = 1.0 / std::sqrt(M_PI * max_freq * mu_si * sigma_si);
            target_dx_skin = skin_depth / 2.0;
        }

        double target_dx = std::min(target_dx_wave, target_dx_skin);
        int auto_segments = static_cast<int>(std::ceil(layer_width / target_dx));

        layer.num_of_segments = std::clamp(auto_segments, 4, 1000);

        // Manual overrides (backward compatibility with old JSON format)
        if (mat_json.contains("num_of_segments")) {
            layer.num_of_segments = mat_json["num_of_segments"].get<size_t>();
        }
        if (mat_json.contains("order")) {
            layer.order = mat_json["order"].get<size_t>();
        }

        // Compute n_skin_depths for all ranks (used for CFL relaxation)
        if (sigma_si > 0.0 && skin_depth > 0.0) {
            layer.n_skin_depths = layer_width / skin_depth;
        }

        if (Mpi::WorldRank() == 0) {
            std::cout << "\n[SGBC Layer Auto-Mesh]" << std::endl;
            std::cout << "  Width                : " << layer_width * 1000.0 << " mm" << std::endl;
            std::cout << "  Pulse Max Freq       : " << max_freq / 1e9 << " GHz" << std::endl;
            std::cout << "  Wavelength           : " << wavelength * 1000.0 << " mm" << std::endl;
            if (sigma_si > 0.0) {
                std::cout << "  Skin Depth           : " << skin_depth * 1000.0 << " mm" << std::endl;
                std::cout << "  Loss Tangent (Pe)    : " << peclet_number << std::endl;
            }
            std::cout << "  Generated Mesh       : " << layer.num_of_segments << " segments (Order " << layer.order << ")" << std::endl;

            if (layer.num_of_segments == 1000) {
                std::cout << "  [WARNING] Max segment limit reached!" << std::endl;
            }
            if (sigma_si > 0.0 && skin_depth > 0.0) {
                double n_skin_depths = layer.n_skin_depths;
                if (n_skin_depths > 7.0) {
                    double transmission_dB = -20.0 * n_skin_depths * std::log10(std::exp(1.0));
                    std::ostringstream oss;
                    oss << "Layer (sigma=" << sigma_si << " S/m, width="
                        << layer_width * 1000.0 << " mm): "
                        << std::fixed << std::setprecision(1)
                        << n_skin_depths << std::defaultfloat
                        << " skin depths (" << std::fixed << std::setprecision(0)
                        << transmission_dB << std::defaultfloat
                        << " dB). Consider PEC.";
                    sgbc_notices.push_back(oss.str());
                }
            }
            std::cout << std::endl;
        }

        return layer;
    };

    for (int b = 0; b < case_data["model"]["boundaries"].size(); b++) {
        if (case_data["model"]["boundaries"][b].contains("type") && 
            case_data["model"]["boundaries"][b]["type"] == "SGBC") { 

            const auto& bdr_json = case_data["model"]["boundaries"][b];

            SGBCProperties props;

            for (auto t = 0; t < bdr_json["tags"].size(); t++) {
                props.geom_tags.emplace_back(bdr_json["tags"][t]);
            }

            // Support both single "material" and multi-layer "layers" formats
            if (bdr_json.contains("layers")) {
                for (const auto& layer_json : bdr_json["layers"]) {
                    props.layers.push_back(parseSGBCLayer(layer_json));
                }
            } else if (bdr_json.contains("material")) {
                props.layers.push_back(parseSGBCLayer(bdr_json["material"]));
            } else {
                throw std::runtime_error("SGBC boundary must define either 'material' or 'layers'. Verify .json parameters.");
            }

            // Parse sgbc_boundaries from boundary level or from material level (backward compat)
            SGBCBoundaryInfo left;
            SGBCBoundaryInfo right;
            auto parseBoundaries = [&](const nlohmann::json& src) {
                if (src.contains("sgbc_boundaries")) {
                    if (src["sgbc_boundaries"].contains("left")) {
                        left.isOn = true;
                        left.bdrCond = assignBoundaryType(src["sgbc_boundaries"]["left"]);
                    }
                    if (src["sgbc_boundaries"].contains("right")) {
                        right.isOn = true;
                        right.bdrCond = assignBoundaryType(src["sgbc_boundaries"]["right"]);
                    }
                }
            };
            parseBoundaries(bdr_json);
            if (!left.isOn && !right.isOn && bdr_json.contains("material")) {
                parseBoundaries(bdr_json["material"]);
            }
            props.sgbc_bdr_info = std::make_pair(left, right);

            if (Mpi::WorldRank() == 0) {
                std::cout << "[SGBC] Total: " << props.layers.size() << " layer(s), "
                          << props.totalSegments() << " segments, "
                          << props.totalWidth() * 1000.0 << " mm width" << std::endl;
            }

            sgbc_props.emplace_back(props);
        }
    }
    
    res.setSGBCProperties(sgbc_props);

    if (Mpi::WorldRank() == 0 && !sgbc_notices.empty()) {
        std::cout << "\n========================================================" << std::endl;
        std::cout << "  SGBC SETUP NOTICES" << std::endl;
        std::cout << "========================================================" << std::endl;
        for (const auto& notice : sgbc_notices) {
            std::cout << "  * " << notice << std::endl;
        }
        std::cout << "========================================================\n" << std::endl;
    }

    return res;
}

json parseJSONfile(const std::string& case_name)
{
	std::ifstream test_file(case_name);
	return json::parse(test_file);
}

void validateCaseNamingStyle(const std::string& case_json_path, const json& case_data)
{
    const std::filesystem::path json_path(case_json_path);
    const std::string json_extension = json_path.extension().string();
    if (json_extension != ".json") {
        throw std::runtime_error(
            "Input case file must be a .json file. Received: " + case_json_path);
    }

    const std::filesystem::path case_dir = json_path.parent_path();
    if (case_dir.empty()) {
        throw std::runtime_error(
            "Invalid case path style for " + case_json_path +
            ". Expected a folder named as the case containing <casename>.json and <casename>.msh.");
    }

    if (!case_data.contains("model") || !case_data["model"].contains("filename") ||
        !case_data["model"]["filename"].is_string()) {
        throw std::runtime_error(
            "Invalid JSON style in " + case_json_path +
            ". model.filename is required and must be a string.");
    }

    const std::string folder_name = case_dir.filename().string();
    const std::string json_name = json_path.stem().string();
    const std::string expected_mesh_name = json_name + ".msh";
    const std::string model_filename_raw = case_data["model"]["filename"].get<std::string>();
    const std::filesystem::path model_filename_path(model_filename_raw);
    const std::string model_filename_name = model_filename_path.filename().string();

    std::ostringstream err;
    err << "Invalid case naming style for " << case_json_path << ". "
        << "Folder name, .json name, .msh name and model.filename must match.\n"
        << "Expected: folder='" << json_name << "', json='" << json_name
        << ".json', msh='" << expected_mesh_name << "', model.filename='"
        << expected_mesh_name << "'.\n"
        << "Found: folder='" << folder_name << "', json='" << json_path.filename().string()
        << "', model.filename='" << model_filename_raw << "'.\n"
        << "Fix the case names so all four values use the same <casename>.";

    if (folder_name != json_name) {
        throw std::runtime_error(err.str());
    }

    if (model_filename_name != expected_mesh_name) {
        throw std::runtime_error(err.str());
    }

    if (model_filename_path.has_parent_path()) {
        throw std::runtime_error(err.str());
    }

    const std::filesystem::path expected_mesh_path = case_dir / expected_mesh_name;
    if (!std::filesystem::exists(expected_mesh_path)) {
        throw std::runtime_error(err.str());
    }
}

maxwell::Solver buildSolverJson(const std::string& case_name, const bool isTest)
{
	auto case_data = parseJSONfile(case_name);
    validateCaseNamingStyle(case_name, case_data);

	return buildSolver(case_data, case_name, isTest);
}

void postProcessInformation(const json& case_data, maxwell::Model& model, maxwell::SolverOptions& solverOpts) 
{
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		mfem::Array<int> tfsf_tags;
		if (case_data["sources"][s]["type"] == "planewave" || case_data["sources"][s]["type"] == "dipole") {
			for (auto t{ 0 }; t < case_data["sources"][s]["tags"].size(); t++) {
				tfsf_tags.Append(case_data["sources"][s]["tags"][t]);
			}
			auto tfsf_atts_present_in_partition_marker{ model.getMarker(maxwell::BdrCond::TotalFieldIn, true) };
			tfsf_atts_present_in_partition_marker.SetSize(model.getConstMesh().bdr_attributes.Max());
			tfsf_atts_present_in_partition_marker = 0;
			for (auto t = 0; t < tfsf_tags.Size(); t++){
				for (auto b = 0; b < model.getConstMesh().GetNBE(); b++){	
					if (model.getMesh().GetBdrAttribute(b) == tfsf_tags[t]){
						tfsf_atts_present_in_partition_marker[model.getMesh().GetBdrAttribute(b) - 1] = 1;
					}
				}
			}
			if (tfsf_atts_present_in_partition_marker.Sum() != 0){
				model.getTotalFieldScatteredFieldToMarker().insert(std::make_pair(maxwell::BdrCond::TotalFieldIn, tfsf_atts_present_in_partition_marker));
			}
		}
	}

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		mfem::Array<int> sgbc_tags = getSGBCTags(case_data);
		if (case_data["model"]["boundaries"][b]["type"] == "SGBC") {
			auto sgbc_atts_present_in_partition_marker{ model.getMarker(maxwell::BdrCond::SGBC, true) };
			sgbc_atts_present_in_partition_marker.SetSize(model.getConstMesh().bdr_attributes.Max());
			sgbc_atts_present_in_partition_marker = 0;
			for (auto t = 0; t < sgbc_tags.Size(); t++){
				for (auto bn = 0; bn < model.getConstMesh().GetNBE(); bn++){	
					if (model.getMesh().GetBdrAttribute(bn) == sgbc_tags[t]){
						sgbc_atts_present_in_partition_marker[model.getMesh().GetBdrAttribute(bn) - 1] = 1;
					}
				}
			}
			if (sgbc_atts_present_in_partition_marker.Sum() != 0){
				model.getSGBCToMarker().insert(std::make_pair(maxwell::BdrCond::SGBC, sgbc_atts_present_in_partition_marker));
			}
		}
	}

	if (model.getBoundaryToMarker().find(BdrCond::SMA) != model.getBoundaryToMarker().end() && solverOpts.evolution.alpha == 0.0 && solverOpts.evolution.op == EvolutionOperatorType::Hesthaven) {
		throw std::runtime_error("Centered SMA with Hesthaven Evolution Operator not supported yet.");
	}

	model.getMesh().SetAttributes();
}


std::string getRunModeTag()
{
    std::string backend;
    if (mfem::Device::Allows(mfem::Backend::CUDA)){
        backend = "cuda-";
        backend.append(std::to_string(Mpi::WorldSize()));
        return backend;
    }
    else{
        if (Mpi::WorldSize() == 1){
            return "single-core";
        }
        else{
            backend = "mpi-";
            backend.append(std::to_string(Mpi::WorldSize()));
            return backend;
        }
    }
}

void prepareExportDirectories(Model& model)
{
	MPI_Comm comm = model.getMesh().GetComm();
    int comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	if (world_rank == 0) {
		std::filesystem::path simExpPath("Exports/" + getRunModeTag() + "/" + model.meshName_ + "/SimulationStats/");
		
		if (std::filesystem::exists(simExpPath)) {
			std::filesystem::remove_all(simExpPath);
		}

		std::filesystem::create_directories(simExpPath);
	}

    MPI_Barrier(comm);
	
}

maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest)
{
	
	maxwell::SolverOptions solverOpts{ buildSolverOptions(case_data) };
	maxwell::Probes probes{ buildProbes(case_data) };
	// Model (mesh) must be built before sources so that auto-mean computation
	// can inspect the TFSF surface geometry to determine the correct delay.
	maxwell::Model model{ buildModel(case_data, case_path, isTest) };
	maxwell::Sources sources{ buildSources(case_data, &model.getConstMesh()) };

	postProcessInformation(case_data, model, solverOpts);
	prepareExportDirectories(model);

	return maxwell::Solver(model, probes, sources, solverOpts);
}

}
