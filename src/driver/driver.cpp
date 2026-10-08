#include "driver.h"
#include "solver/Checkpoint.h"
#include "string"
#include "components/PMLProperties.h"
#include "components/DebyeProperties.h"
#include "components/LorentzProperties.h"

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
#include <optional>
#include <unordered_set>
#include <iterator>
#include <exception>
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
            if (source.contains("type") && source["type"] == "coaxial_port") {
                const double spread = source.value("spread", 1.0);
                if (spread > 0.0 && spread < min_spread) {
                    min_spread = spread;
                    found_gaussian = true;
                }
                continue;
            }
            if (!source.contains("magnitude")) continue;
            const auto& mag = source["magnitude"];
            // delta_gap stores magnitude as a number, not a planewave object.
            if (!mag.is_object()) {
                continue;
            }
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

static std::vector<std::pair<int,int>> buildAllTwoElementInteriorBoundaryPairs(mfem::Mesh& mesh)
{
    std::vector<std::pair<int,int>> res;
    res.reserve(mesh.GetNBE());
    for (int b = 0; b < mesh.GetNBE(); b++) {
        if (auto* f_trans = mesh.GetInternalBdrFaceTransformations(b)) {
            res.emplace_back(f_trans->Elem1No, f_trans->Elem2No);
        }
    }
    return res;
}

static std::vector<std::pair<int,int>>
gatherWeightedConstraintPairs(mfem::Mesh& mesh,
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

// Count how many constrained interior-boundary pairs are currently assigned
// to each rank.
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

// Pick the rank with the lowest weighted load, optionally excluding one rank.
static int leastLoadedRank(const std::vector<long long>& load, int excluded = -1)
{
    int best = -1;
    for (int r = 0; r < static_cast<int>(load.size()); ++r) {
        if (r == excluded) continue;
        if (best < 0 || load[r] < load[best]) best = r;
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

// Verify TFSF TF/SF co-location after DSU partitioning.
// Checks every tagged TFSF face pair directly, then every connected component
// in the TFSF-only element graph (covers multi-TF / multi-SF corner cases).
static void verifyTFSFPartitionColocation(
    mfem::Mesh& mesh,
    const int* partitioning,
    const mfem::Array<int>& tfsf_tags)
{
    if (tfsf_tags.Size() == 0 || Mpi::WorldRank() != 0) return;

    auto tfsf_pairs = buildTwoElementPairsByTagToSort(mesh, tfsf_tags);
    if (tfsf_pairs.empty()) return;

    int split_pairs = 0;
    for (const auto& pr : tfsf_pairs) {
        if (partitioning[pr.first] != partitioning[pr.second]) {
            ++split_pairs;
            if (split_pairs <= 5) {
                std::cout << "[Partition] TFSF SPLIT pair elems "
                          << pr.first << "/" << pr.second << " ranks "
                          << partitioning[pr.first] << "/"
                          << partitioning[pr.second] << "\n";
            }
        }
    }

    const int NE = mesh.GetNE();
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

    std::vector<int> tfsf_degree(NE, 0);
    for (const auto& pr : tfsf_pairs) {
        unite(pr.first, pr.second);
        if (0 <= pr.first && pr.first < NE) ++tfsf_degree[pr.first];
        if (0 <= pr.second && pr.second < NE) ++tfsf_degree[pr.second];
    }

    std::unordered_map<int, std::vector<int>> tfsf_components;
    for (const auto& pr : tfsf_pairs) {
        for (int e : {pr.first, pr.second}) {
            if (0 <= e && e < NE) {
                tfsf_components[find(e)].push_back(e);
            }
        }
    }
    for (auto& [root, elems] : tfsf_components) {
        std::sort(elems.begin(), elems.end());
        elems.erase(std::unique(elems.begin(), elems.end()), elems.end());
    }

    int split_components = 0;
    int max_component = 0;
    int multi_face_elems = 0;
    int max_elem_degree = 0;
    for (int e = 0; e < NE; ++e) {
        if (tfsf_degree[e] > 0) {
            max_elem_degree = std::max(max_elem_degree, tfsf_degree[e]);
            if (tfsf_degree[e] > 1) ++multi_face_elems;
        }
    }

    for (const auto& [root, elems] : tfsf_components) {
        max_component = std::max(max_component, static_cast<int>(elems.size()));
        const int r0 = partitioning[elems.front()];
        for (int e : elems) {
            if (partitioning[e] != r0) {
                ++split_components;
                break;
            }
        }
    }

    std::cout << "[Partition] TFSF verify: pairs=" << tfsf_pairs.size()
              << " split_pairs=" << split_pairs
              << " tfsf_components=" << tfsf_components.size()
              << " split_components=" << split_components
              << " max_component_elems=" << max_component
              << " max_elem_tfsf_faces=" << max_elem_degree
              << " multi_face_elems=" << multi_face_elems;
    if (split_pairs == 0 && split_components == 0) {
        std::cout << " OK";
    } else {
        std::cout << " FAIL";
    }
    std::cout << std::endl;
}

void applyPairwiseConstraintsPartitioning(mfem::Mesh& mesh,
                                          int* partitioning,
                                          const mfem::Array<int>& tfsf_tags,
                                          const mfem::Array<int>& sgbc_tags)
{
    const int NE = mesh.GetNE();
    const int P  = Mpi::WorldSize();
    if (NE == 0 || P <= 1) return;

    auto pairs = buildAllTwoElementInteriorBoundaryPairs(mesh);
    if (pairs.empty()) return;

    // --- Phase 0: Build weighted load model ---
    auto weighted_pairs = gatherWeightedConstraintPairs(mesh, tfsf_tags, sgbc_tags);
    auto weight = buildElementWeights(NE, weighted_pairs);
    auto load   = computeWeightedLoad(P, NE, partitioning, weight);

    // --- Phase 1: Transitive grouping via DSU (union-find) ---
    // Any two-sided interior boundary face is forced rank-local. At corners,
    // one element may participate in multiple such faces, so we use a DSU to
    // merge elements transitively through shared interior-boundary faces.
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
                  << " tags, SGBC=" << sgbc_tags.Size()
                  << " tags, weighted=" << weighted_pairs.size() << " pairs)\n";
        std::cout << "[Partition] Weighted load: min=" << min_load
                  << " max=" << max_load << " ideal=" << std::fixed
                  << std::setprecision(1) << ideal
                  << " imbalance=" << imbalance << "%\n";
        std::cout << "[Partition] Boundary pairs per rank:";
        for (int r = 0; r < P; ++r) std::cout << " R" << r << "=" << bdr_cnt[r];
        std::cout << std::endl;
    }

    verifyTFSFPartitionColocation(mesh, partitioning, tfsf_tags);
}

// METIS-only: assign each constrained-pair component to one rank (root element's).
static void fixSplitConstraintPairs(
    int* partitioning,
    const std::vector<std::pair<int,int>>& pairs)
{
    if (pairs.empty()) return;

    int max_elem = 0;
    for (const auto& pr : pairs) {
        max_elem = std::max({max_elem, pr.first, pr.second});
    }
    const int NE = max_elem + 1;

    std::vector<int> parent(NE);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if (a != b) parent[b] = a;
    };

    for (const auto& pr : pairs) {
        unite(pr.first, pr.second);
    }

    std::unordered_map<int, std::vector<int>> components;
    for (const auto& pr : pairs) {
        for (int e : {pr.first, pr.second}) {
            components[find(e)].push_back(e);
        }
    }
    for (auto& [root, elems] : components) {
        std::sort(elems.begin(), elems.end());
        elems.erase(std::unique(elems.begin(), elems.end()), elems.end());
        const int rank = partitioning[root];
        for (int e : elems) {
            partitioning[e] = rank;
        }
    }

    int n_split = 0;
    for (const auto& pr : pairs) {
        if (partitioning[pr.first] != partitioning[pr.second]) ++n_split;
    }
    if (Mpi::WorldRank() == 0) {
        std::cout << "[Partition] Co-located " << components.size()
                  << " TFSF/SGBC constraint components after METIS";
        if (n_split > 0) {
            std::cout << " (WARNING: " << n_split << " pairs still split)";
        }
        std::cout << "\n";
    }
}

// Assign every element in each constrained pair component to a fixed rank.
static void pinConstraintPairsToRank(
    int* partitioning,
    const std::vector<std::pair<int,int>>& pairs,
    int target_rank)
{
    if (pairs.empty()) return;

    int max_elem = 0;
    for (const auto& pr : pairs) {
        max_elem = std::max({max_elem, pr.first, pr.second});
    }
    const int NE = max_elem + 1;

    std::vector<int> parent(NE);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if (a != b) parent[b] = a;
    };

    for (const auto& pr : pairs) {
        unite(pr.first, pr.second);
    }

    std::unordered_map<int, std::vector<int>> components;
    for (const auto& pr : pairs) {
        for (int e : {pr.first, pr.second}) {
            components[find(e)].push_back(e);
        }
    }
    for (auto& [root, elems] : components) {
        std::sort(elems.begin(), elems.end());
        elems.erase(std::unique(elems.begin(), elems.end()), elems.end());
        for (int e : elems) {
            partitioning[e] = target_rank;
        }
    }
}

static std::unordered_set<int> elementsInPairs(
    const std::vector<std::pair<int,int>>& pairs)
{
    std::unordered_set<int> elems;
    for (const auto& pr : pairs) {
        elems.insert(pr.first);
        elems.insert(pr.second);
    }
    return elems;
}

// Move non-pinned elements off pin_rank to neighbor ranks (or least-loaded).
static int rebalanceExcessFromRank(
    mfem::Mesh& mesh,
    int* partitioning,
    const std::unordered_set<int>& pinned_elements,
    int pin_rank,
    const std::vector<long long>& weight)
{
    const int NE = mesh.GetNE();
    const int P = Mpi::WorldSize();
    if (P <= 1 || pin_rank < 0 || pin_rank >= P) return 0;

    auto load = computeWeightedLoad(P, NE, partitioning, weight);
    const long long total = std::accumulate(load.begin(), load.end(), 0LL);
    const long long ideal = (total + P - 1) / P;

    mesh.ElementToElementTable();
    const mfem::Table& e2e = mesh.ElementToElementTable();

    auto neighborScore = [&](int elem) {
        int score = 0;
        const int n = e2e.RowSize(elem);
        const int* neigh = e2e.GetRow(elem);
        for (int i = 0; i < n; ++i) {
            if (partitioning[neigh[i]] != pin_rank) ++score;
        }
        return score;
    };

    auto bestNeighborTarget = [&](int elem) -> int {
        std::vector<int> neigh_count(P, 0);
        const int n = e2e.RowSize(elem);
        const int* neigh = e2e.GetRow(elem);
        for (int i = 0; i < n; ++i) {
            const int r = partitioning[neigh[i]];
            if (r != pin_rank && 0 <= r && r < P) {
                neigh_count[r]++;
            }
        }
        int best_r = -1;
        int best_cnt = 0;
        for (int r = 0; r < P; ++r) {
            if (r == pin_rank) continue;
            if (neigh_count[r] > best_cnt) {
                best_cnt = neigh_count[r];
                best_r = r;
            }
        }
        if (best_r >= 0 && best_cnt > 0) return best_r;
        return leastLoadedRank(load, pin_rank);
    };

    int moved = 0;
    while (load[pin_rank] > ideal) {
        int best_elem = -1;
        int best_score = -1;
        for (int e = 0; e < NE; ++e) {
            if (partitioning[e] != pin_rank) continue;
            if (pinned_elements.count(e)) continue;
            const int score = neighborScore(e);
            if (score > best_score) {
                best_score = score;
                best_elem = e;
            }
        }
        if (best_elem < 0) break;

        const int target = bestNeighborTarget(best_elem);
        if (target < 0 || target == pin_rank) break;

        load[pin_rank] -= weight[best_elem];
        load[target] += weight[best_elem];
        partitioning[best_elem] = target;
        ++moved;
    }
    return moved;
}

static void applyMetisPartitioningWithTFSFPinRank0(
    mfem::Mesh& mesh,
    int* partitioning,
    const mfem::Array<int>& tfsf_tags,
    const mfem::Array<int>& sgbc_tags)
{
    const int P = Mpi::WorldSize();
    constexpr int tfsf_pin_rank = 0;

    std::vector<std::pair<int,int>> tfsf_pairs;
    if (tfsf_tags.Size() > 0) {
        tfsf_pairs = buildTwoElementPairsByTagToSort(mesh, tfsf_tags);
    }
    std::vector<std::pair<int,int>> sgbc_pairs;
    if (sgbc_tags.Size() > 0) {
        sgbc_pairs = buildTwoElementPairsByTagToSort(mesh, sgbc_tags);
    }

    const auto weighted_pairs = gatherWeightedConstraintPairs(mesh, tfsf_tags, sgbc_tags);
    const auto weight = buildElementWeights(mesh.GetNE(), weighted_pairs);

    if (!tfsf_pairs.empty()) {
        pinConstraintPairsToRank(partitioning, tfsf_pairs, tfsf_pin_rank);
    }
    if (!sgbc_pairs.empty()) {
        fixSplitConstraintPairs(partitioning, sgbc_pairs);
    }

    const auto pinned_tfsf = elementsInPairs(tfsf_pairs);
    const auto load_after_pin = computeWeightedLoad(P, mesh.GetNE(), partitioning, weight);
    const int moved = rebalanceExcessFromRank(
        mesh, partitioning, pinned_tfsf, tfsf_pin_rank, weight);
    auto load_after = computeWeightedLoad(P, mesh.GetNE(), partitioning, weight);

    if (Mpi::WorldRank() == 0) {
        const long long total = std::accumulate(load_after.begin(), load_after.end(), 0LL);
        const long long ideal = (total + P - 1) / P;
        long long max_load = *std::max_element(load_after.begin(), load_after.end());
        long long min_load = *std::min_element(load_after.begin(), load_after.end());
        double imbalance = (ideal > 0)
            ? (max_load - static_cast<double>(ideal)) / ideal * 100.0
            : 0.0;

        std::cout << "[Partition] TFSF pinned to rank " << tfsf_pin_rank
                  << " (" << pinned_tfsf.size() << " elements)\n";
        std::cout << "[Partition] Rebalanced " << moved
                  << " non-TFSF elements off rank " << tfsf_pin_rank
                  << " (load R" << tfsf_pin_rank << ": "
                  << load_after_pin[tfsf_pin_rank] << " -> "
                  << load_after[tfsf_pin_rank] << ", ideal="
                  << ideal << ")\n";
        std::cout << "[Partition] Weighted load: min=" << min_load
                  << " max=" << max_load << " ideal=" << ideal
                  << " imbalance=" << std::fixed << std::setprecision(1)
                  << imbalance << "%\n";
        std::cout << "[Partition] Elements per rank:";
        for (int r = 0; r < P; ++r) {
            int cnt = 0;
            for (int e = 0; e < mesh.GetNE(); ++e) {
                if (partitioning[e] == r) ++cnt;
            }
            std::cout << " R" << r << "=" << cnt;
        }
        std::cout << std::endl;
    }

    verifyTFSFPartitionColocation(mesh, partitioning, tfsf_tags);
}

static void applyMetisPartitioningWithPairFix(
    mfem::Mesh& mesh,
    int* partitioning,
    const mfem::Array<int>& tfsf_tags,
    const mfem::Array<int>& sgbc_tags)
{
    auto pairs = gatherWeightedConstraintPairs(mesh, tfsf_tags, sgbc_tags);
    if (!pairs.empty()) {
        fixSplitConstraintPairs(partitioning, pairs);
    }
    if (Mpi::WorldRank() == 0) {
        std::cout << "[Partition] METIS-only (DSU rebalance disabled)\n";
    }
    verifyTFSFPartitionColocation(mesh, partitioning, tfsf_tags);
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
	if (input.size() != 3) {
		throw std::runtime_error("Expected a 3-vector.");
	}
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

// Sum of tagged boundary-element lengths. Parallel meshes split faces across ranks.
static double deltaGapCurveLength(const mfem::Mesh& mesh, const json& tags)
{
	double length = 0.0;
	const int dim = mesh.Dimension();
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		if (!bdrElemHasTag(mesh, be, tags)) {
			continue;
		}
		mfem::Array<int> verts;
		mesh.GetBdrElementVertices(be, verts);
		for (int i = 1; i < verts.Size(); ++i) {
			const double* a = mesh.GetVertex(verts[i - 1]);
			const double* b = mesh.GetVertex(verts[i]);
			double s = 0.0;
			for (int d = 0; d < dim; ++d) {
				const double diff = a[d] - b[d];
				s += diff * diff;
			}
			length += std::sqrt(s);
		}
	}
	double global_length = 0.0;
	MPI_Allreduce(&length, &global_length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return global_length;
}

// Parallel-plate sizes of a rectangular 3D gap. h is the vertex span along the
// electric polarization. w is the span in the face, perpendicular to it.
// Spans are global min/max, so a face split across ranks is not summed twice.
struct DeltaGapPlate {
	double separation = 0.0;
	double width = 0.0;
};

static DeltaGapPlate deltaGapPlateSpans(
	const mfem::Mesh& mesh, const json& tags, const mfem::Vector& polarization)
{
	double e[3] = {polarization[0], polarization[1], polarization[2]};
	const double en = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
	e[0] /= en;
	e[1] /= en;
	e[2] /= en;

	double nsum[3] = {0.0, 0.0, 0.0};
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		if (!bdrElemHasTag(mesh, be, tags)) {
			continue;
		}
		// MFEM caches the face transformation on a non-const Mesh.
		mfem::ElementTransformation* T =
			const_cast<mfem::Mesh&>(mesh).GetBdrElementTransformation(be);
		mfem::IntegrationPoint ip;
		ip.Set2(0.0, 0.0);
		T->SetIntPoint(&ip);
		mfem::Vector nor(3);
		mfem::CalcOrtho(T->Jacobian(), nor);
		const double nn = nor.Norml2();
		if (!(nn > 0.0)) {
			continue;
		}
		double n[3] = {nor(0) / nn, nor(1) / nn, nor(2) / nn};
		for (int d = 0; d < 3; ++d) {
			if (std::abs(n[d]) > 1e-12) {
				if (n[d] < 0.0) {
					n[0] = -n[0];
					n[1] = -n[1];
					n[2] = -n[2];
				}
				break;
			}
		}
		nsum[0] = n[0];
		nsum[1] = n[1];
		nsum[2] = n[2];
		break;
	}
	double nall[3] = {0.0, 0.0, 0.0};
	MPI_Allreduce(nsum, nall, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	const double nrm = std::sqrt(nall[0] * nall[0] + nall[1] * nall[1] + nall[2] * nall[2]);
	if (!(nrm > 0.0)) {
		throw std::runtime_error("delta_gap face has a zero normal.");
	}
	nall[0] /= nrm;
	nall[1] /= nrm;
	nall[2] /= nrm;

	const double t[3] = {
		nall[1] * e[2] - nall[2] * e[1],
		nall[2] * e[0] - nall[0] * e[2],
		nall[0] * e[1] - nall[1] * e[0]
	};
	const double tn = std::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
	if (!(tn > 1e-8)) {
		throw std::runtime_error("delta_gap polarization is not tangent to the gap face.");
	}

	const double absent = std::numeric_limits<double>::max();
	double hmin = absent;
	double hmax = -absent;
	double wmin = absent;
	double wmax = -absent;
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		if (!bdrElemHasTag(mesh, be, tags)) {
			continue;
		}
		mfem::Array<int> verts;
		mesh.GetBdrElementVertices(be, verts);
		for (int i = 0; i < verts.Size(); ++i) {
			const double* p = mesh.GetVertex(verts[i]);
			double along = 0.0;
			double across = 0.0;
			for (int d = 0; d < 3; ++d) {
				along += p[d] * e[d];
				across += p[d] * t[d];
			}
			hmin = std::min(hmin, along);
			hmax = std::max(hmax, along);
			wmin = std::min(wmin, across);
			wmax = std::max(wmax, across);
		}
	}
	const double local[4] = {hmin, -hmax, wmin, -wmax};
	double glob[4] = {0.0, 0.0, 0.0, 0.0};
	MPI_Allreduce(local, glob, 4, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	DeltaGapPlate plate;
	plate.separation = -glob[1] - glob[0];
	plate.width = -glob[3] - glob[2];
	if (!(plate.separation > 0.0) || !(plate.width > 0.0)
		|| !std::isfinite(plate.separation) || !std::isfinite(plate.width)) {
		throw std::runtime_error("delta_gap rectangle has no positive separation and width.");
	}
	return plate;
}

struct CoaxialPortGeometry {
	mfem::Vector center;
	mfem::Vector axis;
	double inner_radius = 0.0;
	double outer_radius = 0.0;
	int total_field_volume = -1;
	int scattered_field_volume = -1;
};

static bool jsonListHas(const json& tags, int attr)
{
	if (!tags.is_array()) {
		return false;
	}
	for (const auto& t : tags) {
		if (t.get<int>() == attr) {
			return true;
		}
	}
	return false;
}

static mfem::Vector elementCentroid(const mfem::Mesh& mesh, int el)
{
	mfem::Array<int> verts;
	mesh.GetElementVertices(el, verts);
	const int dim = mesh.Dimension();
	mfem::Vector c(3);
	c = 0.0;
	if (verts.Size() == 0) {
		return c;
	}
	for (int v = 0; v < verts.Size(); ++v) {
		const double* p = mesh.GetVertex(verts[v]);
		for (int d = 0; d < dim && d < 3; ++d) {
			c[d] += p[d];
		}
	}
	c /= static_cast<double>(verts.Size());
	return c;
}

static double radiusToAxis(const double* p, const double* center, const double* axis)
{
	double rel[3];
	double along = 0.0;
	for (int d = 0; d < 3; ++d) {
		rel[d] = p[d] - center[d];
		along += rel[d] * axis[d];
	}
	double perp2 = 0.0;
	for (int d = 0; d < 3; ++d) {
		const double radial = rel[d] - along * axis[d];
		perp2 += radial * radial;
	}
	return std::sqrt(std::max(0.0, perp2));
}

static json collectSmaTags(const json& case_data)
{
	json tags = json::array();
	if (!case_data.contains("model") || !case_data["model"].contains("boundaries")) {
		return tags;
	}
	for (const auto& boundary : case_data["model"]["boundaries"]) {
		if (!boundary.contains("type") || boundary["type"] != "SMA") {
			continue;
		}
		for (const auto& tag : boundary["tags"]) {
			tags.push_back(tag);
		}
	}
	return tags;
}

static json collectPmlVolumeTags(const json& case_data)
{
	json tags = json::array();
	if (!case_data.contains("model") || !case_data["model"].contains("materials")) {
		return tags;
	}
	for (const auto& material : case_data["model"]["materials"]) {
		if (!material.contains("type") || material["type"] != "PML" || !material.contains("tags")) {
			continue;
		}
		for (const auto& tag : material["tags"]) {
			tags.push_back(tag);
		}
	}
	return tags;
}

static std::vector<int> allgatherInts(const std::vector<int>& local)
{
	int nprocs = 1;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	const int rank_count = static_cast<int>(local.size());
	std::vector<int> counts(static_cast<std::size_t>(nprocs), 0);
	MPI_Allgather(&rank_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	std::vector<int> displs(static_cast<std::size_t>(nprocs), 0);
	int total = 0;
	for (int i = 0; i < nprocs; ++i) {
		displs[static_cast<std::size_t>(i)] = total;
		total += counts[static_cast<std::size_t>(i)];
	}
	std::vector<int> gathered(static_cast<std::size_t>(std::max(total, 1)));
	const int dummy = 0;
	MPI_Allgatherv(
		local.empty() ? &dummy : local.data(), rank_count, MPI_INT,
		gathered.data(), counts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
	gathered.resize(static_cast<std::size_t>(std::max(total, 0)));
	return gathered;
}

static std::vector<double> allgatherDoubles(const std::vector<double>& local)
{
	int nprocs = 1;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	const int rank_count = static_cast<int>(local.size());
	std::vector<int> counts(static_cast<std::size_t>(nprocs), 0);
	MPI_Allgather(&rank_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	std::vector<int> displs(static_cast<std::size_t>(nprocs), 0);
	int total = 0;
	for (int i = 0; i < nprocs; ++i) {
		displs[static_cast<std::size_t>(i)] = total;
		total += counts[static_cast<std::size_t>(i)];
	}
	std::vector<double> gathered(static_cast<std::size_t>(std::max(total, 1)));
	const double dummy = 0.0;
	MPI_Allgatherv(
		local.empty() ? &dummy : local.data(), rank_count, MPI_DOUBLE,
		gathered.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
	gathered.resize(static_cast<std::size_t>(std::max(total, 0)));
	return gathered;
}

static double boundaryFaceDistance(const mfem::Mesh& mesh, int be, const double center[3])
{
	mfem::Array<int> verts;
	mesh.GetBdrElementVertices(be, verts);
	if (verts.Size() == 0) {
		return std::numeric_limits<double>::infinity();
	}
	double face[3] = {0.0, 0.0, 0.0};
	const int dim = std::min(mesh.Dimension(), 3);
	for (int v = 0; v < verts.Size(); ++v) {
		const double* p = mesh.GetVertex(verts[v]);
		for (int d = 0; d < dim; ++d) {
			face[d] += p[d];
		}
	}
	const double inv = 1.0 / static_cast<double>(verts.Size());
	double dist2 = 0.0;
	for (int d = 0; d < 3; ++d) {
		const double diff = face[d] * inv - center[d];
		dist2 += diff * diff;
	}
	return std::sqrt(dist2);
}

static int agreeVolumeAttribute(int local, const char* name)
{
	const int has = local >= 0 ? 1 : 0;
	int any = 0;
	MPI_Allreduce(&has, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (any == 0) {
		throw std::runtime_error(std::string("coaxial_port could not identify the ") + name + " volume.");
	}
	const int send_min = has ? local : std::numeric_limits<int>::max();
	const int send_max = has ? local : -1;
	int global_min = 0;
	int global_max = 0;
	MPI_Allreduce(&send_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&send_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (global_min != global_max) {
		throw std::runtime_error(std::string("coaxial_port ") + name + " volume disagrees across ranks.");
	}
	return global_min;
}

static CoaxialPortGeometry fitCoaxialPort(
	mfem::Mesh& mesh,
	const json& load_tags,
	const json& live_tags,
	const json& outer_tags,
	const json& sma_tags,
	const json& pml_tags)
{
	std::unordered_set<int> load_volumes;
	std::vector<int> load_elem_pairs;
	std::vector<int> boundary_volumes;
	int local_n_boundary = 0;
	int local_n_interior = 0;
	int local_boundary_not_sma = 0;
	double load_min[3] = {
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max()};
	double load_max[3] = {
		-std::numeric_limits<double>::max(),
		-std::numeric_limits<double>::max(),
		-std::numeric_limits<double>::max()};
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		if (!jsonListHas(load_tags, mesh.GetBdrAttribute(be))) {
			continue;
		}
		mfem::Array<int> verts;
		mesh.GetBdrElementVertices(be, verts);
		for (int v = 0; v < verts.Size(); ++v) {
			const double* p = mesh.GetVertex(verts[v]);
			for (int d = 0; d < mesh.Dimension() && d < 3; ++d) {
				load_min[d] = std::min(load_min[d], p[d]);
				load_max[d] = std::max(load_max[d], p[d]);
			}
		}
		auto* tr = mesh.GetInternalBdrFaceTransformations(be);
		const bool interior = tr != nullptr && tr->Elem2No >= 0;
		if (!interior) {
			++local_n_boundary;
			if (!jsonListHas(sma_tags, mesh.GetBdrAttribute(be))) {
				local_boundary_not_sma = 1;
			}
			if (tr != nullptr && tr->Elem1No >= 0) {
				boundary_volumes.push_back(mesh.GetAttribute(tr->Elem1No));
			}
			else {
				int el = -1;
				int info = 0;
				mesh.GetBdrElementAdjacentElement(be, el, info);
				if (el >= 0) {
					boundary_volumes.push_back(mesh.GetAttribute(el));
				}
			}
			continue;
		}
		++local_n_interior;
		const int elem1_volume = mesh.GetAttribute(tr->Elem1No);
		const int elem2_volume = mesh.GetAttribute(tr->Elem2No);
		load_volumes.insert(elem1_volume);
		load_volumes.insert(elem2_volume);
		load_elem_pairs.push_back(elem1_volume);
		load_elem_pairs.push_back(elem2_volume);
	}
	int global_n_boundary = 0;
	int global_n_interior = 0;
	int global_boundary_not_sma = 0;
	MPI_Allreduce(&local_n_boundary, &global_n_boundary, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_n_interior, &global_n_interior, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_boundary_not_sma, &global_boundary_not_sma, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (global_n_boundary > 0 && global_n_interior > 0) {
		throw std::runtime_error("coaxial_port load tags must be either all interior or all boundary faces.");
	}
	if (global_n_boundary == 0 && global_n_interior == 0) {
		throw std::runtime_error("coaxial_port load tags match no faces.");
	}
	const bool boundary_port = global_n_boundary > 0;
	if (boundary_port && global_boundary_not_sma != 0) {
		throw std::runtime_error("coaxial_port boundary load must also be an SMA boundary.");
	}

	const int local_too_many = load_volumes.size() > 2 ? 1 : 0;
	int global_too_many = 0;
	MPI_Allreduce(&local_too_many, &global_too_many, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (global_too_many != 0) {
		throw std::runtime_error("coaxial_port load face touches more than two volumes.");
	}

	double center[3] = {0.0, 0.0, 0.0};
	for (int d = 0; d < 3; ++d) {
		double global_min = 0.0;
		double global_max = 0.0;
		MPI_Allreduce(&load_min[d], &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&load_max[d], &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		center[d] = 0.5 * (global_min + global_max);
	}

	std::unordered_set<int> pml_volumes;
	for (const auto& tag : pml_tags) {
		pml_volumes.insert(tag.get<int>());
	}

	std::vector<int> local_sma;
	std::vector<int> local_edges;
	std::vector<int> local_interface_volume;
	std::vector<double> local_interface_distance;
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		auto* tr = mesh.GetInternalBdrFaceTransformations(be);
		if (jsonListHas(sma_tags, mesh.GetBdrAttribute(be))) {
			if (tr != nullptr && tr->Elem1No >= 0) {
				local_sma.push_back(mesh.GetAttribute(tr->Elem1No));
				if (tr->Elem2No >= 0) {
					local_sma.push_back(mesh.GetAttribute(tr->Elem2No));
				}
			}
			else {
				int el = -1;
				int info = 0;
				mesh.GetBdrElementAdjacentElement(be, el, info);
				if (el >= 0) {
					local_sma.push_back(mesh.GetAttribute(el));
				}
			}
		}
		if (tr == nullptr || tr->Elem2No < 0 || jsonListHas(load_tags, mesh.GetBdrAttribute(be))) {
			continue;
		}
		const int side_a = mesh.GetAttribute(tr->Elem1No);
		const int side_b = mesh.GetAttribute(tr->Elem2No);
		local_edges.push_back(side_a);
		local_edges.push_back(side_b);
		const bool a_pml = pml_volumes.count(side_a) != 0;
		const bool b_pml = pml_volumes.count(side_b) != 0;
		if (a_pml == b_pml) {
			continue;
		}
		local_interface_volume.push_back(a_pml ? side_b : side_a);
		local_interface_distance.push_back(boundaryFaceDistance(mesh, be, center));
	}

	const std::vector<int> load_ids = allgatherInts(
		std::vector<int>(load_volumes.begin(), load_volumes.end()));
	const std::vector<int> load_pairs = allgatherInts(load_elem_pairs);
	const std::vector<int> boundary_ids = allgatherInts(boundary_volumes);
	const std::vector<int> sma_ids = allgatherInts(local_sma);
	const std::vector<int> edges = allgatherInts(local_edges);
	const std::vector<int> interface_volumes = allgatherInts(local_interface_volume);
	const std::vector<double> interface_distances = allgatherDoubles(local_interface_distance);

	std::unordered_set<int> global_load(load_ids.begin(), load_ids.end());
	std::unordered_set<int> global_sma(sma_ids.begin(), sma_ids.end());
	std::unordered_map<int, std::vector<int>> adjacent;
	for (std::size_t i = 0; i + 1 < edges.size(); i += 2) {
		adjacent[edges[i]].push_back(edges[i + 1]);
		adjacent[edges[i + 1]].push_back(edges[i]);
	}
	std::unordered_map<int, double> interface_distance;
	for (std::size_t i = 0; i < interface_volumes.size() && i < interface_distances.size(); ++i) {
		const int volume = interface_volumes[i];
		const double distance = interface_distances[i];
		const auto found = interface_distance.find(volume);
		if (found == interface_distance.end() || distance < found->second) {
			interface_distance[volume] = distance;
		}
	}

	int local_tf = -1;
	int local_sf = -1;
	if (boundary_port) {
		std::unordered_set<int> boundary_set(boundary_ids.begin(), boundary_ids.end());
		if (boundary_set.size() != 1) {
			throw std::runtime_error(
				"coaxial_port boundary load must lie on exactly one volume.");
		}
		local_tf = *boundary_set.begin();
	}
	else if (global_load.size() == 2) {
		std::vector<int> sides(global_load.begin(), global_load.end());
		std::sort(sides.begin(), sides.end());
		const bool first_sma = global_sma.count(sides[0]) != 0;
		const bool second_sma = global_sma.count(sides[1]) != 0;
		if (first_sma != second_sma) {
			local_sf = first_sma ? sides[0] : sides[1];
			local_tf = first_sma ? sides[1] : sides[0];
		}
		else if (!first_sma) {
			const auto pmlDistance = [&](int start) {
				std::unordered_set<int> seen;
				std::vector<int> stack;
				stack.push_back(start);
				seen.insert(start);
				double best = std::numeric_limits<double>::infinity();
				bool reached = false;
				if (pml_volumes.count(start) != 0) {
					reached = true;
					best = 0.0;
				}
				for (std::size_t n = 0; n < stack.size(); ++n) {
					const int volume = stack[n];
					const auto iface = interface_distance.find(volume);
					if (iface != interface_distance.end()) {
						reached = true;
						best = std::min(best, iface->second);
					}
					if (pml_volumes.count(volume) != 0) {
						reached = true;
					}
					for (const int next : adjacent[volume]) {
						if (seen.insert(next).second) {
							stack.push_back(next);
						}
					}
				}
				return reached ? best : std::numeric_limits<double>::infinity();
			};
			const double first_distance = pmlDistance(sides[0]);
			const double second_distance = pmlDistance(sides[1]);
			const bool first_reached = std::isfinite(first_distance);
			const bool second_reached = std::isfinite(second_distance);
			if (first_reached != second_reached) {
				local_sf = first_reached ? sides[0] : sides[1];
				local_tf = first_reached ? sides[1] : sides[0];
			}
			else if (first_reached && second_reached) {
				const double scale = std::max(1.0, std::max(first_distance, second_distance));
				if (!(std::abs(first_distance - second_distance) > 1.0e-8 * scale)) {
					throw std::runtime_error(
						"coaxial_port could not tell the scattered-field volume from the total-field volume. "
						"One side of the load must meet an SMA boundary.");
				}
				const bool first_nearer = first_distance < second_distance;
				local_sf = first_nearer ? sides[0] : sides[1];
				local_tf = first_nearer ? sides[1] : sides[0];
			}
			else {
				throw std::runtime_error(
					"coaxial_port could not tell the scattered-field volume from the total-field volume. "
					"One side of the load must meet an SMA boundary.");
			}
		}
		else {
			if (load_pairs.size() < 2 || (load_pairs.size() % 2) != 0) {
				throw std::runtime_error(
					"coaxial_port could not read the load face Elem1/Elem2 order.");
			}
			const int elem1_volume = load_pairs[0];
			const int elem2_volume = load_pairs[1];
			for (std::size_t i = 0; i + 1 < load_pairs.size(); i += 2) {
				if (load_pairs[i] != elem1_volume || load_pairs[i + 1] != elem2_volume) {
					throw std::runtime_error(
						"coaxial_port load face Elem1/Elem2 order is not the same on every face.");
				}
			}
			local_sf = elem1_volume;
			local_tf = elem2_volume;
		}
	}
	const int tf_volume = agreeVolumeAttribute(local_tf, "total-field");
	const int sf_volume = boundary_port ? -1 : agreeVolumeAttribute(local_sf, "scattered-field");

	double sum_tf[3] = {0.0, 0.0, 0.0};
	double sum_sf[3] = {0.0, 0.0, 0.0};
	int n_tf = 0;
	int n_sf = 0;
	for (int el = 0; el < mesh.GetNE(); ++el) {
		const int attr = mesh.GetAttribute(el);
		if (attr != tf_volume && attr != sf_volume) {
			continue;
		}
		const mfem::Vector c = elementCentroid(mesh, el);
		double* sum = attr == tf_volume ? sum_tf : sum_sf;
		int& n = attr == tf_volume ? n_tf : n_sf;
		for (int d = 0; d < 3; ++d) {
			sum[d] += c[d];
		}
		++n;
	}
	double global_sum_tf[3];
	double global_sum_sf[3];
	int global_n_tf = 0;
	int global_n_sf = 0;
	MPI_Allreduce(sum_tf, global_sum_tf, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(sum_sf, global_sum_sf, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&n_tf, &global_n_tf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&n_sf, &global_n_sf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (global_n_tf == 0 || (!boundary_port && global_n_sf == 0)) {
		throw std::runtime_error("coaxial_port total-field or scattered-field volume has no elements.");
	}
	double axis[3];
	double axis_norm = 0.0;
	for (int d = 0; d < 3; ++d) {
		if (boundary_port) {
			axis[d] = global_sum_tf[d] / global_n_tf - center[d];
		}
		else {
			axis[d] = global_sum_tf[d] / global_n_tf - global_sum_sf[d] / global_n_sf;
		}
		axis_norm += axis[d] * axis[d];
	}
	axis_norm = std::sqrt(axis_norm);
	if (!(axis_norm > 0.0)) {
		throw std::runtime_error(
			boundary_port
				? "coaxial_port boundary load and its volume share a center."
				: "coaxial_port total-field and scattered-field volumes share a centroid.");
	}
	for (int d = 0; d < 3; ++d) {
		axis[d] /= axis_norm;
	}

	auto meanCylinderRadius = [&](const json& tags) {
		double sum = 0.0;
		int count = 0;
		for (const auto& tag_json : tags) {
			const int tag = tag_json.get<int>();
			std::unordered_set<int> seen;
			double min_r = std::numeric_limits<double>::max();
			double max_r = 0.0;
			double tag_sum = 0.0;
			int tag_count = 0;
			for (int be = 0; be < mesh.GetNBE(); ++be) {
				if (mesh.GetBdrAttribute(be) != tag) {
					continue;
				}
				mfem::Array<int> verts;
				mesh.GetBdrElementVertices(be, verts);
				for (int v = 0; v < verts.Size(); ++v) {
					if (!seen.insert(verts[v]).second) {
						continue;
					}
					const double r = radiusToAxis(mesh.GetVertex(verts[v]), center, axis);
					min_r = std::min(min_r, r);
					max_r = std::max(max_r, r);
					tag_sum += r;
					++tag_count;
				}
			}
			double global_min = 0.0;
			double global_max = 0.0;
			double global_sum = 0.0;
			int global_count = 0;
			const double send_min = tag_count > 0 ? min_r : std::numeric_limits<double>::max();
			MPI_Allreduce(&send_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			MPI_Allreduce(&max_r, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&tag_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&tag_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (global_count == 0 || !(global_max > 0.0)) {
				continue;
			}
			if ((global_max - global_min) / global_max > 0.15) {
				continue;
			}
			sum += global_sum;
			count += global_count;
		}
		if (count == 0) {
			return 0.0;
		}
		return sum / static_cast<double>(count);
	};

	const double inner = meanCylinderRadius(live_tags);
	const double outer = meanCylinderRadius(outer_tags);
	if (!(inner > 0.0) || !(outer > inner)) {
		throw std::runtime_error(
			"coaxial_port could not measure inner and outer radii from the cylindrical conductor tags.");
	}

	CoaxialPortGeometry geom;
	geom.center.SetSize(3);
	geom.axis.SetSize(3);
	for (int d = 0; d < 3; ++d) {
		geom.center[d] = center[d];
		geom.axis[d] = axis[d];
	}
	geom.inner_radius = inner;
	geom.outer_radius = outer;
	geom.total_field_volume = tf_volume;
	geom.scattered_field_volume = sf_volume;
	return geom;
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
	const double gaussMean,
	const double amplitude_peak = 1.0,
	const double peak_radius = 1.0)
{
	DerivGaussDipole dip(length, gaussianSpread, gaussMean, amplitude_peak, peak_radius);
	return std::make_unique<TotalField>(dip);
}

Sources buildSources(const json& case_data, const mfem::Mesh* mesh)
{
	Sources res;
	int delta_gap_count = 0;
	int coaxial_count = 0;
	bool closed_tfsf = false;
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
			if (coaxial_count > 0) {
				throw std::runtime_error(
					"coaxial_port is the TFSF source. Do not combine it with planewave or dipole.");
			}
			closed_tfsf = true;
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
			if (coaxial_count > 0) {
				throw std::runtime_error(
					"coaxial_port is the TFSF source. Do not combine it with planewave or dipole.");
			}
			closed_tfsf = true;
			const auto& mag = case_data["sources"][s]["magnitude"];
			if (mag.contains("amplitude")) {
				throw std::runtime_error(
					"Dipole magnitude.amplitude was renamed to magnitude.amplitude_peak "
					"(desired max |E| at peak_radius).");
			}
			double length = mag["length"].get<double>();
			double spread = mag["spread"].get<double>();
			double amplitude_peak = mag.value("amplitude_peak", 1.0);
			if (amplitude_peak <= 0.0) {
				throw std::runtime_error(
					"Dipole magnitude.amplitude_peak must be > 0 (got " +
					std::to_string(amplitude_peak) + ").");
			}

			// Determine mean: use explicit value if provided, otherwise auto-compute.
			// For the retarded-time Gaussian, mean_auto ensures the field is
			// AUTO_DELAY_N_SIGMA sigma before peak at t=0 on the TFSF surface.
			double mean;
			double peak_radius = 1.0;
			if (mag.contains("peak_radius")) {
				peak_radius = mag["peak_radius"].get<double>();
				if (peak_radius <= 0.0) {
					throw std::runtime_error(
						"Dipole magnitude.peak_radius must be > 0.");
				}
			}
			if (mag.contains("mean")) {
				mean = mag["mean"].get<double>();
			} else if (mesh && case_data["sources"][s].contains("tags")) {
				double min_r = minRadiusOnTFSFSurface(*mesh, case_data["sources"][s]["tags"]);
				// mean_auto = N*sigma*sqrt(2) - r_min/c
				// (r_min/c is the retarded-time advance the TFSF surface already provides)
				double mean_auto = AUTO_DELAY_N_SIGMA * spread * std::sqrt(2.0)
				                   - min_r / physicalConstants::speedOfLight;
				mean = std::max(spread * std::sqrt(2.0), mean_auto);
				if (!mag.contains("peak_radius")) {
					peak_radius = min_r;
				}
				std::cout << "[Source " << s << " (dipole)] Auto-computed mean = "
				          << mean << " (min_radius=" << min_r
				          << "), amplitude_peak=" << amplitude_peak
				          << " at peak_radius=" << peak_radius << "\n";
			} else {
				mean = AUTO_DELAY_N_SIGMA * spread * std::sqrt(2.0);
			}

			res.add(buildDerivGaussDipole(
				length, spread, mean, amplitude_peak, peak_radius));
		}
		else if (case_data["sources"][s]["type"] == "delta_gap") {
			if (delta_gap_count > 0) {
				throw std::runtime_error("Only one delta_gap source is supported.");
			}
			++delta_gap_count;
			if (mesh == nullptr) {
				throw std::runtime_error("delta_gap requires a mesh.");
			}
			if (mesh->Dimension() != 2 && mesh->Dimension() != 3) {
				throw std::runtime_error("delta_gap requires a 2D or 3D mesh.");
			}
			const auto& src = case_data["sources"][s];
			if (!src.contains("tags") || src["tags"].empty()) {
				throw std::runtime_error("delta_gap requires tags.");
			}
			if (src.contains("normal_sign")) {
				throw std::runtime_error(
					"delta_gap no longer accepts normal_sign. Set polarization, the electric-field direction.");
			}
			if (!src.contains("polarization") || !src["polarization"].is_array()
				|| src["polarization"].size() != 3) {
				throw std::runtime_error("delta_gap polarization must be a 3-vector.");
			}
			mfem::Vector polarization = assemble3DVector(src["polarization"]);
			for (int d = 0; d < polarization.Size(); ++d) {
				if (!std::isfinite(polarization[d])) {
					throw std::runtime_error("delta_gap polarization must be finite.");
				}
			}
			if (!(polarization.Norml2() > 0.0)) {
				throw std::runtime_error("delta_gap polarization must be nonzero.");
			}
			double magnitude = 1.0;
			if (src.contains("magnitude")) {
				if (!src["magnitude"].is_number()) {
					throw std::runtime_error("delta_gap magnitude must be a number.");
				}
				magnitude = src["magnitude"].get<double>();
			}
			if (!(magnitude > 0.0) || !std::isfinite(magnitude)) {
				throw std::runtime_error("delta_gap magnitude must be > 0.");
			}
			double db_cut = -20.0;
			if (src.contains("db_cut")) {
				if (!src["db_cut"].is_number()) {
					throw std::runtime_error("delta_gap db_cut must be a number.");
				}
				db_cut = src["db_cut"].get<double>();
			}
			bool derivative = false;
			if (src.contains("signal")) {
				if (!src["signal"].is_string()) {
					throw std::runtime_error("delta_gap signal must be a string.");
				}
				const std::string signal = src["signal"].get<std::string>();
				if (signal == "gaussian") {
					derivative = false;
				}
				else if (signal == "gaussian_derivative") {
					derivative = true;
				}
				else {
					throw std::runtime_error(
						"delta_gap signal must be \"gaussian\" or \"gaussian_derivative\".");
				}
			}
			if (src.contains("f_max")) {
				throw std::runtime_error(
					"delta_gap no longer accepts f_max. Set spread, the Gaussian width in normalized time.");
			}
			const double curve_length = deltaGapCurveLength(*mesh, src["tags"]);
			if (!(curve_length > 0.0)) {
				throw std::runtime_error("delta_gap tags match no curve on the mesh.");
			}
			double spread = 0.0;
			if (src.contains("spread")) {
				if (src.contains("db_cut")) {
					throw std::runtime_error(
						"delta_gap spread sets the pulse width directly. Omit db_cut.");
				}
				if (!src["spread"].is_number()) {
					throw std::runtime_error("delta_gap spread must be a number.");
				}
				spread = src["spread"].get<double>();
			}
			else {
				spread = gaussianSpreadForDbCut(magnitude * curve_length / 10.0, db_cut);
			}
			if (!(spread > 0.0) || !std::isfinite(spread)) {
				throw std::runtime_error("delta_gap spread must be > 0.");
			}
			const double t0 = AUTO_DELAY_N_SIGMA * spread * std::sqrt(2.0);
			const bool rectangular = mesh->Dimension() == 3;
			DeltaGapPlate plate;
			if (rectangular) {
				plate = deltaGapPlateSpans(*mesh, src["tags"], polarization);
			}
			int rank = 0;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0) {
				std::cout << "[delta_gap] L=" << curve_length;
				if (rectangular) {
					// Wide-plate TEM values, eta = 1. Diagnostics, not source scales.
					const double impedance = plate.separation / plate.width;
					const double capacitance = plate.width / plate.separation;
					std::cout << " h=" << plate.separation
					          << " w=" << plate.width
					          << " Z=" << impedance
					          << " C'=" << capacitance
					          << " w/h=" << (plate.width / plate.separation);
				}
				std::cout << " magnitude=" << magnitude
				          << " signal=" << (derivative ? "gaussian_derivative" : "gaussian")
				          << " spread=" << spread
				          << " t0=" << t0 << "\n";
			}
			res.add(std::make_unique<DeltaGapSource>(magnitude, spread, t0, polarization, derivative));
		}
		else if (case_data["sources"][s]["type"] == "coaxial_port") {
			if (coaxial_count > 0) {
				throw std::runtime_error("Only one coaxial_port source is supported.");
			}
			++coaxial_count;
			if (closed_tfsf) {
				throw std::runtime_error(
					"coaxial_port is the TFSF source. Do not combine it with planewave or dipole.");
			}
			if (mesh == nullptr) {
				throw std::runtime_error("coaxial_port requires a mesh.");
			}
			if (mesh->Dimension() != 3) {
				throw std::runtime_error("coaxial_port requires a 3D mesh.");
			}
			const auto& src = case_data["sources"][s];
			if (!src.contains("tags") || !src["tags"].is_object()) {
				throw std::runtime_error("coaxial_port tags must name outer, live, and load.");
			}
			for (const char* key : {"outer", "live", "load"}) {
				if (!src["tags"].contains(key) || !src["tags"][key].is_array() || src["tags"][key].empty()) {
					throw std::runtime_error(
						std::string("coaxial_port tags.") + key + " must be a non-empty array.");
				}
			}
			double magnitude = 1.0;
			if (src.contains("magnitude")) {
				if (!src["magnitude"].is_number()) {
					throw std::runtime_error("coaxial_port magnitude must be a number.");
				}
				magnitude = src["magnitude"].get<double>();
			}
			if (!(magnitude > 0.0) || !std::isfinite(magnitude)) {
				throw std::runtime_error("coaxial_port magnitude must be > 0.");
			}
			double spread = 1.0;
			if (src.contains("spread")) {
				if (!src["spread"].is_number()) {
					throw std::runtime_error("coaxial_port spread must be a number.");
				}
				spread = src["spread"].get<double>();
			}
			if (!(spread > 0.0) || !std::isfinite(spread)) {
				throw std::runtime_error("coaxial_port spread must be > 0.");
			}
			bool derivative = false;
			if (src.contains("signal")) {
				if (!src["signal"].is_string()) {
					throw std::runtime_error("coaxial_port signal must be a string.");
				}
				const std::string signal = src["signal"].get<std::string>();
				if (signal == "gaussian") {
					derivative = false;
				}
				else if (signal == "gaussian_derivative") {
					derivative = true;
				}
				else {
					throw std::runtime_error(
						"coaxial_port signal must be \"gaussian\" or \"gaussian_derivative\".");
				}
			}
			const CoaxialPortGeometry geom = fitCoaxialPort(
				const_cast<mfem::Mesh&>(*mesh), src["tags"]["load"], src["tags"]["live"], src["tags"]["outer"],
				collectSmaTags(case_data), collectPmlVolumeTags(case_data));
			const double t0 = AUTO_DELAY_N_SIGMA * spread * std::sqrt(2.0);
			int rank = 0;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0) {
				const auto previous = std::cout.precision();
				std::cout << std::setprecision(6);
				std::cout << "[coaxial_port] a=" << geom.inner_radius
				          << " b=" << geom.outer_radius
				          << " Z=" << std::log(geom.outer_radius / geom.inner_radius) / (2.0 * std::acos(-1.0))
				          << " tf_volume=" << geom.total_field_volume
				          << " sf_volume=" << geom.scattered_field_volume
				          << " axis=(" << geom.axis[0] << ", " << geom.axis[1] << ", " << geom.axis[2] << ")"
				          << " magnitude=" << magnitude
				          << " spread=" << spread
				          << " t0=" << t0 << "\n";
				std::cout << std::setprecision(static_cast<int>(previous));
			}
			res.add(std::make_unique<TotalField>(CoaxialMode(
				magnitude, spread, t0, derivative,
				geom.center, geom.axis,
				geom.inner_radius, geom.outer_radius,
				geom.total_field_volume, geom.scattered_field_volume)));
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

		if (case_data["solver_options"].contains("checkpoint_percent")) {
			const double percent = case_data["solver_options"]["checkpoint_percent"].get<double>();
			if (!(percent >= 0.0 && percent <= 100.0)) {
				throw std::runtime_error(
					"solver_options.checkpoint_percent must be between 0 and 100.");
			}
			res.setCheckpointPercent(percent);
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
            const auto& exp_json = case_data["probes"]["exporter"];
            if (exp_json.contains("saves")) {
                throw std::runtime_error(
                    "probes.exporter: \"saves\" was replaced by \"save_every\" "
                    "(solver-time interval). Example: \"save_every\": 0.5 with "
                    "final_time 20 writes t=0,0.5,...,20.");
            }
            if (exp_json.contains("save_every") && exp_json.contains("steps")) {
                throw std::runtime_error(
                    "probes.exporter: specify either \"save_every\" or \"steps\", not both.");
            }
            ExporterProbe exporter_probe;
            if (exp_json.contains("name")) {
                exporter_probe.name = exp_json["name"];
            } else {
                exporter_probe.name = case_data["model"]["filename"];
            }
            if (exp_json.contains("save_every")) {
                exporter_probe.save_every = exp_json["save_every"].get<double>();
                if (!(exporter_probe.save_every > 0.0)) {
                    throw std::runtime_error(
                        "probes.exporter.save_every must be > 0.");
                }
            } else {
                exporter_probe.visSteps = calculate_interval(exp_json);
            }
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

static bool isSGBCBoundaryType(const std::string& boundary_type);

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
    else if (isSGBCBoundaryType(bdr_cond)) {
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

GeomTagToMaterialInfo assembleAttributeToMaterial(
	const json& case_data, const mfem::Mesh& mesh)
{
	GeomTagToMaterialInfo res{};

	checkIfThrows(case_data.contains("model"), "JSON data does not include 'model'.");
	checkIfThrows(case_data["model"].contains("materials"), "JSON data does not include 'materials'.");

	std::unordered_map<int, std::string> material_tag_kind;
	auto claimMaterialTag = [&](int tag, const std::string& kind) {
		const auto it = material_tag_kind.find(tag);
		if (it == material_tag_kind.end()) {
			material_tag_kind.emplace(tag, kind);
			return;
		}
		const bool dispersive =
			kind == "debye" || kind == "lorentz" ||
			it->second == "debye" || it->second == "lorentz";
		const bool pml = kind == "pml" || it->second == "pml";
		if (pml && dispersive) {
			throw std::runtime_error(kDispersiveOnPmlNotAllowed);
		}
		if (dispersive) {
			throw std::runtime_error(
				"Material tag " + std::to_string(tag) +
				" cannot carry Debye or Lorentz together with another material assignment.");
		}
	};

	for (auto m = 0; m < case_data["model"]["materials"].size(); m++) {
		const auto& mat_json = case_data["model"]["materials"][m];
		if (!mat_json.contains("tags")) {
			continue;
		}
		std::string kind = "material";
		if (mat_json.contains("type")) {
			const std::string type = mat_json["type"].get<std::string>();
			if (type == "vacuum") {
				kind = "vacuum";
				if (mat_json.contains("debye") || mat_json.contains("lorentz")) {
					throw std::runtime_error(
						"Vacuum material must not define debye or lorentz.");
				}
			} else if (type == "PML") {
				kind = "pml";
				if (mat_json.contains("debye") || mat_json.contains("lorentz")) {
					throw std::runtime_error(kDispersiveOnPmlNotAllowed);
				}
			}
		} else if (mat_json.contains("debye") && mat_json.contains("lorentz")) {
			throw std::runtime_error(
				"A material cannot define both debye and lorentz.");
		} else if (mat_json.contains("debye")) {
			kind = "debye";
		} else if (mat_json.contains("lorentz")) {
			kind = "lorentz";
		}
		for (auto t = 0; t < mat_json["tags"].size(); t++) {
			claimMaterialTag(mat_json["tags"][t].get<int>(), kind);
		}
	}

	for (auto m = 0; m < case_data["model"]["materials"].size(); m++) {
		const auto& mat_json = case_data["model"]["materials"][m];

		if (mat_json.contains("type")) {
			const std::string type = mat_json["type"].get<std::string>();
			if (type == "vacuum") {
				if (mat_json.contains("debye") || mat_json.contains("lorentz")) {
					throw std::runtime_error(
						"Vacuum material must not define debye or lorentz.");
				}
				const Material vacuum = buildVacuumMaterial();
				for (auto t = 0; t < mat_json["tags"].size(); t++) {
					res.gt2m.emplace(mat_json["tags"][t], vacuum);
				}
			} else if (type == "PML") {
				PMLProperties props =
					parsePMLMaterialBlock(mat_json, mesh.Dimension());
				const Material vacuum = buildVacuumMaterial();
				for (auto t = 0; t < mat_json["tags"].size(); t++) {
					const GeomTag tag = mat_json["tags"][t];
					props.geom_tags.push_back(tag);
					res.gt2m.emplace(tag, vacuum);
				}
				res.pml_props.push_back(std::move(props));
			} else {
				throw std::runtime_error(
					"Unknown material type '" + type + "'. Supported: vacuum, PML, or legacy eps/mu.");
			}
			continue;
		}

		if (mat_json.contains("debye")) {
			if (mat_json.contains("relative_permittivity")) {
				throw std::runtime_error(
					"Debye material must not define relative_permittivity. "
					"The electric mass uses debye.eps_inf.");
			}
			DebyeProperties pole = parseDebyeObject(mat_json["debye"]);
			double mu{ 1.0 }, sigma{ 0.0 };
			if (mat_json.contains("relative_permeability")) {
				mu = mat_json["relative_permeability"];
			}
			if (mat_json.contains("bulk_conductivity")) {
				sigma = mat_json["bulk_conductivity"].get<double>() * physicalConstants::freeSpaceImpedance_SI;
			}
			for (auto t = 0; t < mat_json["tags"].size(); t++) {
				pole.geom_tag = mat_json["tags"][t].get<int>();
				res.debye.push_back(pole);
				res.gt2m.emplace(pole.geom_tag, Material(pole.eps_inf, mu, sigma));
			}
			continue;
		}

		if (mat_json.contains("lorentz")) {
			if (mat_json.contains("debye")) {
				throw std::runtime_error(
					"A material cannot define both debye and lorentz.");
			}
			if (mat_json.contains("relative_permittivity")) {
				throw std::runtime_error(
					"Lorentz material must not define relative_permittivity. "
					"The electric mass uses lorentz.eps_inf.");
			}
			LorentzProperties pole = parseLorentzObject(mat_json["lorentz"]);
			double mu{ 1.0 }, sigma{ 0.0 };
			if (mat_json.contains("relative_permeability")) {
				mu = mat_json["relative_permeability"];
			}
			if (mat_json.contains("bulk_conductivity")) {
				sigma = mat_json["bulk_conductivity"].get<double>() * physicalConstants::freeSpaceImpedance_SI;
			}
			for (auto t = 0; t < mat_json["tags"].size(); t++) {
				pole.geom_tag = mat_json["tags"][t].get<int>();
				res.lorentz.push_back(pole);
				res.gt2m.emplace(pole.geom_tag, Material(pole.eps_inf, mu, sigma));
			}
			continue;
		}

		for (auto t = 0; t < mat_json["tags"].size(); t++) {
			double eps{ 1.0 }, mu{ 1.0 }, sigma{ 0.0 };
			if (mat_json.contains("relative_permittivity")) {
				eps = mat_json["relative_permittivity"];
			}
			if (mat_json.contains("relative_permeability")) {
				mu = mat_json["relative_permeability"];
			}
			if (mat_json.contains("bulk_conductivity")) {
				sigma = mat_json["bulk_conductivity"].get<double>() * physicalConstants::freeSpaceImpedance_SI;
			}
			res.gt2m.emplace(std::make_pair(mat_json["tags"][t], Material(eps, mu, sigma)));
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
	int failed = 0;
	std::string err;
	if (Mpi::WorldRank() == 0) {
		try {
			fixGmshMesh(mesh_string);
		} catch (const std::exception& ex) {
			failed = 1;
			err = ex.what();
		}
	}
	int global_failed = 0;
	MPI_Allreduce(&failed, &global_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (global_failed) {
		if (!err.empty()) {
			std::cerr << err << std::endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
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
                res.Append(case_data["sources"][s]["tags"][t].get<int>());
            }
        }
        else if (case_data["sources"][s].contains("type") &&
                 case_data["sources"][s]["type"] == "coaxial_port" &&
                 case_data["sources"][s].contains("tags") &&
                 case_data["sources"][s]["tags"].contains("load")) {
            for (int t = 0; t < case_data["sources"][s]["tags"]["load"].size(); ++t) {
                res.Append(case_data["sources"][s]["tags"]["load"][t].get<int>());
            }
        }
    }
    return res;
}

static mfem::Array<int> getDeltaGapTags(const json& case_data)
{
    mfem::Array<int> res;
    if (!case_data.contains("sources")) {
        return res;
    }
    for (int s = 0; s < case_data["sources"].size(); ++s) {
        if (!case_data["sources"][s].contains("type")
            || case_data["sources"][s]["type"] != "delta_gap"
            || !case_data["sources"][s].contains("tags")) {
            continue;
        }
        const auto& tags = case_data["sources"][s]["tags"];
        for (int t = 0; t < tags.size(); ++t) {
            res.Append(tags[t].get<int>());
        }
    }
    return res;
}

static bool isSGBCBoundaryType(const std::string& boundary_type)
{
	return boundary_type == "SGBC";
}

Array<int> getSGBCTags(const json& case_data)
{
    Array<int> res;
    if (!case_data.contains("model") ||
        !case_data["model"].contains("boundaries")) {
        return res;
    }

    for (int b = 0; b < case_data["model"]["boundaries"].size(); ++b) {
        if (case_data["model"]["boundaries"][b].contains("type") &&
            isSGBCBoundaryType(case_data["model"]["boundaries"][b]["type"].get<std::string>())) {
            for (int t = 0; t < case_data["model"]["boundaries"][b]["tags"].size(); ++t) {
                res.Append(case_data["model"]["boundaries"][b]["tags"][t].get<int>());
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

Model buildModel(
    const json& case_data,
    const std::string& case_path,
    const bool isTest,
    const int* partition_override,
    int partition_count)
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

    mfem::Array<int> tfsf_tags = getTFSFTags(case_data);
    mfem::Array<int> sgbc_tags = getSGBCTags(case_data);
    {
        mfem::Array<int> gap_tags = getDeltaGapTags(case_data);
        for (int i = 0; i < gap_tags.Size(); ++i) {
            sgbc_tags.Append(gap_tags[i]);
        }
    }

    const int ne = mesh.GetNE();
    mfem::Array<int> part(ne);
    int partition_bad = 0;
    std::string partition_error;
    if (Mpi::WorldRank() == 0) {
        if (partition_override != nullptr) {
            if (partition_count != ne) {
                partition_bad = 1;
                partition_error = "Checkpoint partition has " + std::to_string(partition_count)
                    + " elements; the mesh has " + std::to_string(ne) + ".";
            } else {
                for (int i = 0; i < ne; ++i) {
                    if (partition_override[i] < 0 || partition_override[i] >= Mpi::WorldSize()) {
                        partition_bad = 1;
                        partition_error = "Checkpoint partition rank is outside this job.";
                        break;
                    }
                    part[i] = partition_override[i];
                }
            }
        } else {
            int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
            const char* use_dsu = std::getenv("DGTD_USE_DSU_PARTITION");
            const char* pin_tfsf = std::getenv("DGTD_TFSF_PIN_RANK0");
            const bool use_tfsf_pin_rank0 = tfsf_tags.Size() > 0 &&
                (!pin_tfsf || pin_tfsf[0] != '0');

            if (use_dsu && use_dsu[0] == '1') {
                applyPairwiseConstraintsPartitioning(mesh, partitioning, tfsf_tags, sgbc_tags);
            } else if (use_tfsf_pin_rank0) {
                applyMetisPartitioningWithTFSFPinRank0(mesh, partitioning, tfsf_tags, sgbc_tags);
            } else {
                applyMetisPartitioningWithPairFix(mesh, partitioning, tfsf_tags, sgbc_tags);
            }
            for (int i = 0; i < ne; ++i) {
                part[i] = partitioning[i];
            }
            delete[] partitioning;
        }
    }
    MPI_Bcast(&partition_bad, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (partition_bad) {
        if (!partition_error.empty()) {
            throw std::runtime_error(partition_error);
        }
        throw std::runtime_error("Checkpoint partition does not match this job.");
    }
    if (Mpi::WorldSize() > 1 && ne > 0) {
        MPI_Bcast(part.GetData(), ne, MPI_INT, 0, MPI_COMM_WORLD);
    }

    Model res(mesh, att_to_material, att_to_bdr_info, ne > 0 ? part.GetData() : nullptr);
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
        if (mat_json.contains("debye") || mat_json.contains("lorentz")) {
            throw std::runtime_error(
                "SGBC layer must not define debye or lorentz.");
        }
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
            isSGBCBoundaryType(case_data["model"]["boundaries"][b]["type"].get<std::string>())) {

            const auto& bdr_json = case_data["model"]["boundaries"][b];

            SGBCProperties props;

            for (auto t = 0; t < bdr_json["tags"].size(); t++) {
                props.geom_tags.emplace_back(bdr_json["tags"][t]);
            }
            props.exporter_probe = bdr_json.value("exporter_probe", false);

            // Optional solver_options.sgbc_cfl replaces the historical hard-coded 0.5
            // in recommended_dt_ = crossing_time * sgbc_cfl * opacity_relax.
            if (case_data.contains("solver_options") &&
                case_data["solver_options"].contains("sgbc_cfl")) {
                const double cfl = case_data["solver_options"]["sgbc_cfl"].get<double>();
                if (!(cfl > 0.0)) {
                    throw std::runtime_error("solver_options.sgbc_cfl must be > 0.");
                }
                props.sgbc_cfl = cfl;
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

    res.setDebyeProperties(att_to_material.debye);
    res.setLorentzProperties(att_to_material.lorentz);
    res.setPMLProperties(att_to_material.pml_props);
    if (res.hasPML() && Mpi::WorldRank() == 0) {
        std::cout << "\n[PML] Parsed " << att_to_material.pml_props.size() << " region(s):" << std::endl;
        for (size_t ri = 0; ri < att_to_material.pml_props.size(); ++ri) {
            const auto& props = att_to_material.pml_props[ri];
            std::cout << "  Region " << ri << ": " << props.geom_tags.size()
                      << " tag(s), stretch_mode="
                      << (props.stretch_mode == PMLStretchMode::Radial ? "radial" : "box")
                      << ", grading_order=" << props.grading_order
                      << ", target_reflection=" << std::scientific << props.target_reflection
                      << std::defaultfloat
                      << ", kappa_max=" << props.kappa_max
                      << ", alpha_max=" << props.alpha_max
                      << ", active_axes:";
            for (Direction d : props.active_axes) {
                std::cout << " " << d;
            }
            std::cout << std::endl;
        }
        // Profiles are initialized in Solver with the case FE order (MPI-safe
        // attribute-based σ eval; serial mesh for global interface / L).
    }

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

static std::string meshCaseNameFromJson(const json& case_data)
{
	const std::string filename = case_data.at("model").at("filename").get<std::string>();
	if (filename.size() >= 5 && filename.compare(filename.size() - 5, 5, ".mesh") == 0) {
		return filename.substr(0, filename.size() - 5);
	}
	if (filename.size() >= 4 && filename.compare(filename.size() - 4, 4, ".msh") == 0) {
		return filename.substr(0, filename.size() - 4);
	}
	throw std::runtime_error("File format for mesh must be name.msh or name.mesh");
}

maxwell::Solver buildSolverJson(const std::string& case_name, const bool isTest, bool restart)
{
	auto case_data = parseJSONfile(case_name);
    validateCaseNamingStyle(case_name, case_data);

	return buildSolver(case_data, case_name, isTest, restart);
}

static int boundaryAttributeSlots(const mfem::Mesh& mesh)
{
	if (mesh.bdr_attributes.Size() == 0) {
		return 0;
	}
	return mesh.bdr_attributes.Max();
}

static bool markLocalBoundaryAttribute(mfem::Array<int>& marker, int attr)
{
	if (attr < 1 || attr > marker.Size()) {
		return false;
	}
	marker[attr - 1] = 1;
	return true;
}

static void throwIfAnyRankFailed(int local_bad, MPI_Comm comm, const char* what)
{
	int global_bad = 0;
	MPI_Allreduce(&local_bad, &global_bad, 1, MPI_INT, MPI_MAX, comm);
	if (global_bad) {
		throw std::runtime_error(
			std::string(what)
			+ " boundary attribute is outside the mesh boundary-attribute range.");
	}
}

void postProcessInformation(const json& case_data, maxwell::Model& model, maxwell::SolverOptions& solverOpts) 
{
	const MPI_Comm comm = model.getMesh().GetComm();

	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		mfem::Array<int> tfsf_tags;
		if (case_data["sources"][s]["type"] == "planewave" || case_data["sources"][s]["type"] == "dipole") {
			for (auto t{ 0 }; t < case_data["sources"][s]["tags"].size(); t++) {
				tfsf_tags.Append(case_data["sources"][s]["tags"][t].get<int>());
			}
		}
		else if (case_data["sources"][s]["type"] == "coaxial_port") {
			const auto& load = case_data["sources"][s]["tags"]["load"];
			for (auto t{ 0 }; t < load.size(); t++) {
				tfsf_tags.Append(load[t].get<int>());
			}
		}
		if (tfsf_tags.Size() == 0) {
			continue;
		}
			const int n_attr = boundaryAttributeSlots(model.getConstMesh());
			if (n_attr <= 0) {
				throw std::runtime_error(
					"TFSF tag is set but the mesh has no boundary attributes.");
			}
			auto tfsf_atts_present_in_partition_marker{ model.getMarker(maxwell::BdrCond::TotalFieldIn, true) };
			tfsf_atts_present_in_partition_marker.SetSize(n_attr);
			tfsf_atts_present_in_partition_marker = 0;
			int tfsf_bad = 0;
			for (auto t = 0; t < tfsf_tags.Size(); t++){
				for (auto b = 0; b < model.getConstMesh().GetNBE(); b++){	
					if (model.getMesh().GetBdrAttribute(b) == tfsf_tags[t]){
						if (!markLocalBoundaryAttribute(
								tfsf_atts_present_in_partition_marker,
								model.getMesh().GetBdrAttribute(b))) {
							tfsf_bad = 1;
						}
					}
				}
			}
			throwIfAnyRankFailed(tfsf_bad, comm, "TFSF");
			const int local_tfsf_marker = tfsf_atts_present_in_partition_marker.Sum();
			int global_tfsf_marker = 0;
			MPI_Allreduce(&local_tfsf_marker, &global_tfsf_marker, 1, MPI_INT, MPI_SUM, comm);
			if (global_tfsf_marker != 0) {
				model.getTotalFieldScatteredFieldToMarker().insert(
					std::make_pair(maxwell::BdrCond::TotalFieldIn, tfsf_atts_present_in_partition_marker));
			}
	}

	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		if (case_data["sources"][s]["type"] != "delta_gap") {
			continue;
		}
		mfem::Array<int> gap_tags;
		for (auto t{ 0 }; t < case_data["sources"][s]["tags"].size(); t++) {
			gap_tags.Append(case_data["sources"][s]["tags"][t].get<int>());
		}
		mfem::Array<int> gap_marker;
		const int n_attr = model.getConstMesh().bdr_attributes.Size() == 0
			? 0 : model.getConstMesh().bdr_attributes.Max();
		gap_marker.SetSize(n_attr);
		gap_marker = 0;
		int gap_bad = 0;
		for (int t = 0; t < gap_tags.Size(); t++) {
			for (int b = 0; b < model.getConstMesh().GetNBE(); b++) {
				if (model.getMesh().GetBdrAttribute(b) == gap_tags[t]) {
					if (!markLocalBoundaryAttribute(
							gap_marker, model.getMesh().GetBdrAttribute(b))) {
						gap_bad = 1;
					}
				}
			}
		}
		throwIfAnyRankFailed(gap_bad, comm, "delta_gap");
		const int local_gap_marker = gap_marker.Sum();
		int global_gap_marker = 0;
		MPI_Allreduce(&local_gap_marker, &global_gap_marker, 1, MPI_INT, MPI_SUM, comm);
		if (global_gap_marker == 0) {
			throw std::runtime_error("delta_gap tags were not found on the mesh.");
		}
		model.setDeltaGapMarker(gap_marker);
	}

    mfem::Array<int> sgbc_tags = getSGBCTags(case_data);
    if (sgbc_tags.Size() != 0) {
        const int n_attr = boundaryAttributeSlots(model.getConstMesh());
        if (n_attr <= 0) {
            throw std::runtime_error(
                "SGBC tag is set but the mesh has no boundary attributes.");
        }
        auto sgbc_atts_present_in_partition_marker{ model.getMarker(maxwell::BdrCond::SGBC, true) };
        sgbc_atts_present_in_partition_marker.SetSize(n_attr);
        sgbc_atts_present_in_partition_marker = 0;
        int sgbc_bad = 0;
        for (auto t = 0; t < sgbc_tags.Size(); t++){
            for (auto bn = 0; bn < model.getConstMesh().GetNBE(); bn++){	
                if (model.getMesh().GetBdrAttribute(bn) == sgbc_tags[t]){
                    if (!markLocalBoundaryAttribute(
                            sgbc_atts_present_in_partition_marker,
                            model.getMesh().GetBdrAttribute(bn))) {
                        sgbc_bad = 1;
                    }
                }
            }
        }
        throwIfAnyRankFailed(sgbc_bad, comm, "SGBC");
        const int local_sgbc_marker = sgbc_atts_present_in_partition_marker.Sum();
        int global_sgbc_marker = 0;
        MPI_Allreduce(&local_sgbc_marker, &global_sgbc_marker, 1, MPI_INT, MPI_SUM, comm);
        if (global_sgbc_marker != 0) {
            model.getSGBCToMarker().insert(
                std::make_pair(maxwell::BdrCond::SGBC, sgbc_atts_present_in_partition_marker));
        }
    }

	if (model.getBoundaryToMarker().find(BdrCond::SMA) != model.getBoundaryToMarker().end() && solverOpts.evolution.alpha == 0.0 && solverOpts.evolution.op == EvolutionOperatorType::Hesthaven) {
		throw std::runtime_error("Centered SMA with Hesthaven Evolution Operator not supported yet.");
	}

	model.getMesh().SetAttributes();
}


void prepareExportDirectories(Model& model, bool preserve_existing)
{
	MPI_Comm comm = model.getMesh().GetComm();
    int comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	if (world_rank == 0) {
		std::filesystem::path simExpPath(maxwell::getSimulationCaseExportPath(model.meshName_) + "/SimulationStats/");
		
		if (std::filesystem::exists(simExpPath) && !preserve_existing) {
			std::filesystem::remove_all(simExpPath);
		}

		std::filesystem::create_directories(simExpPath);
	}

    MPI_Barrier(comm);
	
}

maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest, bool restart)
{
	if (restart && isTest) {
		throw std::runtime_error("--restart is only supported from opensemba_dgtd.");
	}

	maxwell::SolverOptions solverOpts{ buildSolverOptions(case_data) };
	maxwell::Probes probes{ buildProbes(case_data) };

	// If MOR state probes are defined, force export_operator to true
	if (!probes.morStateProbes.empty()) {
		solverOpts.setExportEO(true);
	}

	std::filesystem::path checkpoint_dir;
	std::vector<int> saved_partition;
	std::string mesh_path;
	std::optional<maxwell::CheckpointManifest> checkpoint_manifest;
	if (!isTest) {
		const std::string case_tag = meshCaseNameFromJson(case_data);
		mesh_path = assembleLauncherMeshString(case_data["model"]["filename"].get<std::string>(), case_path);
		const auto latest = maxwell::findLatestCompleteCheckpoint(case_tag);
		const auto foreign = maxwell::findForeignCompleteCheckpoint(case_tag);
		if (!latest && foreign) {
			throw std::runtime_error(maxwell::foreignCheckpointMessage(*foreign));
		}
		if (restart) {
			if (!latest) {
				throw std::runtime_error(
					"No complete checkpoint found under " + maxwell::checkpointRoot(case_tag).string()
					+ ". --restart requires a checkpoint from the same JSON, mesh, and MPI rank count.");
			}
			checkpoint_dir = *latest;
			checkpoint_manifest = maxwell::readManifest(checkpoint_dir / "manifest.json");
			if (checkpoint_manifest->format_version != maxwell::kCheckpointFormatVersion) {
				throw std::runtime_error(
					"Checkpoint format version " + std::to_string(checkpoint_manifest->format_version)
					+ " is not supported (this build reads version "
					+ std::to_string(maxwell::kCheckpointFormatVersion) + ").");
			}
			if (checkpoint_manifest->world_size != mfem::Mpi::WorldSize()) {
				throw std::runtime_error(
					"Checkpoint was written with " + std::to_string(checkpoint_manifest->world_size)
					+ " ranks; this job has " + std::to_string(mfem::Mpi::WorldSize()) + ".");
			}
			if (maxwell::sha256File(case_path) != checkpoint_manifest->json_sha256
				&& !maxwell::sameCaseExceptCheckpointPercent(case_path, checkpoint_dir / "input.json")) {
				throw std::runtime_error(
					"Checkpoint JSON does not match " + case_path
					+ ". Resume requires the same input file. solver_options.checkpoint_percent may change.");
			}
			saved_partition = maxwell::readPartition(checkpoint_dir / "partition.bin");
			solverOpts.resume_from_checkpoint = true;
		} else if (latest) {
			throw std::runtime_error(
				"A complete checkpoint exists at " + latest->string()
				+ ". Resume with --restart, or remove the Checkpoints directory to start from t = 0.");
		}
	}

	// Model (mesh) must be built before sources so that auto-mean computation
	// can inspect the TFSF surface geometry to determine the correct delay.
	const int* partition_override = saved_partition.empty() ? nullptr : saved_partition.data();
	maxwell::Model model{ buildModel(
		case_data,
		case_path,
		isTest,
		partition_override,
		static_cast<int>(saved_partition.size())) };
	if (checkpoint_manifest) {
		if (maxwell::sha256File(mesh_path) != checkpoint_manifest->mesh_sha256) {
			throw std::runtime_error(
				"Checkpoint mesh does not match " + mesh_path + ".");
		}
	}
	maxwell::Sources sources{ buildSources(case_data, &model.getConstMesh()) };

	postProcessInformation(case_data, model, solverOpts);
	prepareExportDirectories(model, restart);

	if (!isTest && solverOpts.checkpoint_percent > 0.0) {
		solverOpts.checkpoint_json_path = case_path;
		solverOpts.checkpoint_mesh_path = mesh_path;
	}
	if (restart) {
		solverOpts.checkpoint_directory = checkpoint_dir.string();
	}
	return maxwell::Solver(model, probes, sources, solverOpts);
}

}
