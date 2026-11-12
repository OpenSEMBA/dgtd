#include "wires/ThinWire.h"

#include "math/Function.h"

#include "math/PhysicalConstants.h"

#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace maxwell {
namespace {

constexpr int kLineType = 1;

struct CurveEdge {
	int tag = 0;
	int a = 0;
	int b = 0;
};

double dot3(const double* a, const double* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void cross3(const double* a, const double* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

double norm3(const double* a)
{
	return std::sqrt(dot3(a, a));
}

bool readGmshCurves(
	const std::string& path,
	std::unordered_map<int, std::array<double, 3>>& nodes,
	std::vector<CurveEdge>& edges)
{
	std::ifstream in(path);
	if (!in) {
		throw std::runtime_error("Cannot open mesh for thin-wire curves: " + path);
	}
	std::string line;
	while (std::getline(in, line)) {
		if (line.find("$Nodes") != std::string::npos) {
			int n = 0;
			in >> n;
			std::getline(in, line);
			for (int i = 0; i < n; ++i) {
				int id;
				std::array<double, 3> x{};
				in >> id >> x[0] >> x[1] >> x[2];
				nodes[id] = x;
			}
		} else if (line.find("$Elements") != std::string::npos) {
			int n = 0;
			in >> n;
			std::getline(in, line);
			for (int i = 0; i < n; ++i) {
				std::string raw;
				if (!std::getline(in, raw)) {
					break;
				}
				std::istringstream iss(raw);
				std::vector<int> tok;
				int value = 0;
				while (iss >> value) {
					tok.push_back(value);
				}
				if (tok.size() < 5) {
					continue;
				}
				const int etype = tok[1];
				if (etype != kLineType) {
					if (etype == 8 || etype == 26 || etype == 27 || etype == 28) {
						throw std::runtime_error(
							"Thin wires accept linear Gmsh line elements (type 1) only.");
					}
					continue;
				}
				const int ntags = tok[2];
				if (ntags < 2 || static_cast<int>(tok.size()) < 3 + ntags + 2) {
					throw std::runtime_error(
						"A thin-wire curve needs an elementary tag on each Gmsh line.");
				}
				CurveEdge edge;
				edge.tag = tok[4];
				edge.a = tok[3 + ntags];
				edge.b = tok[4 + ntags];
				edges.push_back(edge);
			}
		}
	}
	return !edges.empty();
}

WirePolyline chainEdges(
	const std::vector<int>& tags,
	const std::vector<CurveEdge>& edges,
	const std::unordered_map<int, std::array<double, 3>>& nodes)
{
	std::set<int> wanted(tags.begin(), tags.end());
	std::vector<CurveEdge> mine;
	for (const auto& edge : edges) {
		if (wanted.count(edge.tag)) {
			mine.push_back(edge);
		}
	}
	if (mine.empty()) {
		throw std::runtime_error(
			"No Gmsh line carries the elementary tag requested by a wire.");
	}

	std::map<int, std::vector<std::pair<int, int>>> adj;
	for (int e = 0; e < static_cast<int>(mine.size()); ++e) {
		adj[mine[e].a].push_back({e, mine[e].b});
		adj[mine[e].b].push_back({e, mine[e].a});
	}
	for (const auto& item : adj) {
		if (item.second.size() > 2) {
			throw std::runtime_error("Branched thin wires are not supported.");
		}
	}

	std::vector<int> ends;
	for (const auto& item : adj) {
		if (item.second.size() == 1) {
			ends.push_back(item.first);
		}
	}
	const bool loop = ends.empty();
	if (!loop && ends.size() != 2) {
		throw std::runtime_error("A thin wire must be one open curve or one closed loop.");
	}

	std::vector<char> used(mine.size(), 0);
	std::vector<int> node_order;
	int start = loop ? adj.begin()->first : ends[0];
	int prev = -1;
	int current = start;
	node_order.push_back(current);
	for (int step = 0; step < static_cast<int>(mine.size()); ++step) {
		int next_edge = -1;
		int next_node = -1;
		for (const auto& link : adj[current]) {
			if (!used[link.first] && link.second != prev) {
				next_edge = link.first;
				next_node = link.second;
				break;
			}
		}
		if (next_edge < 0) {
			throw std::runtime_error("The thin-wire curve is not a single connected path.");
		}
		used[next_edge] = 1;
		prev = current;
		current = next_node;
		if (!(loop && current == start)) {
			node_order.push_back(current);
		}
		if (loop && current == start) {
			break;
		}
	}
	if (loop) {
		node_order.push_back(start);
	}
	for (char flag : used) {
		if (!flag) {
			throw std::runtime_error("The thin-wire tags do not form one connected curve.");
		}
	}

	WirePolyline poly;
	poly.tags = tags;
	poly.closed = loop;
	poly.points.reserve(node_order.size());
	for (int id : node_order) {
		const auto it = nodes.find(id);
		if (it == nodes.end()) {
			throw std::runtime_error("A thin-wire node id is missing from $Nodes.");
		}
		poly.points.push_back(it->second);
	}
	if (poly.points.size() < 2) {
		throw std::runtime_error("A thin wire needs at least two points.");
	}
	return poly;
}

struct Clip {
	int element = -1;
	double s0 = 0.0;
	double s1 = 1.0;
};

bool clipSegmentToTet(
	const double* a,
	const double* b,
	const double verts[4][3],
	double& s0,
	double& s1)
{
	s0 = 0.0;
	s1 = 1.0;
	const int faces[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
	const int opposite[4] = {0, 1, 2, 3};
	double dir[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
	for (int f = 0; f < 4; ++f) {
		const double* p0 = verts[faces[f][0]];
		const double* p1 = verts[faces[f][1]];
		const double* p2 = verts[faces[f][2]];
		double e1[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
		double e2[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
		double normal[3];
		cross3(e1, e2, normal);
		const double* inside = verts[opposite[f]];
		double toward[3] = {inside[0] - p0[0], inside[1] - p0[1], inside[2] - p0[2]};
		if (dot3(normal, toward) < 0.0) {
			normal[0] = -normal[0];
			normal[1] = -normal[1];
			normal[2] = -normal[2];
		}
		double a_rel[3] = {a[0] - p0[0], a[1] - p0[1], a[2] - p0[2]};
		const double n_a = dot3(normal, a_rel);
		const double n_dir = dot3(normal, dir);
		const double scale = std::max(1.0, norm3(normal));
		if (std::abs(n_dir) <= 1e-14 * scale) {
			if (n_a < -1e-10 * scale) {
				return false;
			}
			continue;
		}
		const double hit = -n_a / n_dir;
		if (n_dir > 0.0) {
			s0 = std::max(s0, hit);
		} else {
			s1 = std::min(s1, hit);
		}
		if (s1 < s0 - 1e-10) {
			return false;
		}
	}
	return s1 >= s0 - 1e-12;
}

double meanEdgeLength(const double verts[4][3])
{
	double sum = 0.0;
	int count = 0;
	for (int i = 0; i < 4; ++i) {
		for (int j = i + 1; j < 4; ++j) {
			double d[3] = {
				verts[j][0] - verts[i][0],
				verts[j][1] - verts[i][1],
				verts[j][2] - verts[i][2]};
			sum += norm3(d);
			++count;
		}
	}
	return sum / static_cast<double>(count);
}

bool markerHits(const mfem::Array<int>& marker, int attribute)
{
	return marker.Size() > 0 && attribute > 0 && attribute <= marker.Size() && marker[attribute - 1] == 1;
}

bool pointOnTriangle(
	const double* p,
	const double* a,
	const double* b,
	const double* c,
	double tol)
{
	double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
	double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
	double ap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
	double normal[3];
	cross3(ab, ac, normal);
	const double area = norm3(normal);
	if (area < 1e-14) {
		return false;
	}
	if (std::abs(dot3(normal, ap)) > tol * area) {
		return false;
	}
	const double d00 = dot3(ab, ab);
	const double d01 = dot3(ab, ac);
	const double d11 = dot3(ac, ac);
	const double d20 = dot3(ap, ab);
	const double d21 = dot3(ap, ac);
	const double denom = d00 * d11 - d01 * d01;
	if (std::abs(denom) < 1e-28) {
		return false;
	}
	const double beta = (d11 * d20 - d01 * d21) / denom;
	const double gamma = (d00 * d21 - d01 * d20) / denom;
	const double alpha = 1.0 - beta - gamma;
	return alpha >= -tol && beta >= -tol && gamma >= -tol;
}

int endCondition(
	const mfem::Mesh& mesh,
	const Model& model,
	const double* point)
{
	const double tol = 1e-8 * std::max(1.0, norm3(point));
	int kind = 0;
	for (int be = 0; be < mesh.GetNBE(); ++be) {
		mfem::Array<int> verts;
		mesh.GetBdrElementVertices(be, verts);
		if (verts.Size() != 3) {
			continue;
		}
		double tri[3][3];
		for (int v = 0; v < 3; ++v) {
			const double* coord = mesh.GetVertex(verts[v]);
			tri[v][0] = coord[0];
			tri[v][1] = coord[1];
			tri[v][2] = coord[2];
		}
		if (!pointOnTriangle(point, tri[0], tri[1], tri[2], tol)) {
			continue;
		}
		const int attr = mesh.GetBdrAttribute(be);
		const auto& boundary = model.getGeomTagToBoundaryCond();
		const auto& interior = model.getGeomTagToIntBoundaryCond();
		auto read = [&](const GeomTagToBoundary& map) {
			const auto it = map.find(attr);
			return it == map.end() ? BdrCond::PEC : it->second;
		};
		BdrCond cond = BdrCond::SMA;
		if (boundary.count(attr)) {
			cond = read(boundary);
		} else if (interior.count(attr)) {
			cond = interior.at(attr);
		} else {
			continue;
		}
		if (cond == BdrCond::PMC || cond == BdrCond::SMA || cond == BdrCond::SGBC) {
			throw std::runtime_error(
				"A thin-wire end on PMC, SMA, or SGBC is not supported.");
		}
		if (cond == BdrCond::PEC) {
			kind = 1;
		}
	}
	return kind;
}

int serialElementToLocal(const Model& model, const mfem::ParMesh& mesh, int serial)
{
	const mfem::Array<int>& part = model.elementPartition();
	if (part.Size() == 0) {
		if (mfem::Mpi::WorldSize() != 1) {
			throw std::runtime_error(
				"Thin wires need the element partition on more than one rank.");
		}
		return serial >= 0 && serial < mesh.GetNE() ? serial : -1;
	}
	if (serial < 0 || serial >= part.Size()) {
		throw std::runtime_error("A thin-wire host element is outside the serial mesh.");
	}
	if (part[serial] != mfem::Mpi::WorldRank()) {
		return -1;
	}
	int local = 0;
	for (int i = 0; i < serial; ++i) {
		if (part[i] == part[serial]) {
			++local;
		}
	}
	if (local >= mesh.GetNE()) {
		throw std::runtime_error("A thin-wire host is missing from the local mesh.");
	}
	return local;
}

mfem::Array<int> uniteMarker(const mfem::Array<int>& marker, MPI_Comm comm, int comm_size)
{
	mfem::Array<int> copy(marker);
	if (comm_size <= 1) {
		return copy;
	}
	const int n = copy.Size();
	int nmax = 0;
	int nmin = 0;
	MPI_Allreduce(&n, &nmax, 1, MPI_INT, MPI_MAX, comm);
	MPI_Allreduce(&n, &nmin, 1, MPI_INT, MPI_MIN, comm);
	if (nmin != nmax) {
		throw std::runtime_error("Thin-wire boundary markers differ in length across ranks.");
	}
	if (n > 0) {
		MPI_Allreduce(MPI_IN_PLACE, copy.GetData(), n, MPI_INT, MPI_MAX, comm);
	}
	return copy;
}

bool segmentMeetsExpandedBox(
	const double* a, const double* b, const double* lo, const double* hi, double pad)
{
	double t0 = 0.0;
	double t1 = 1.0;
	for (int d = 0; d < 3; ++d) {
		const double start = a[d];
		const double dir = b[d] - a[d];
		const double box0 = lo[d] - pad;
		const double box1 = hi[d] + pad;
		if (std::abs(dir) < 1e-15) {
			if (start < box0 || start > box1) {
				return false;
			}
			continue;
		}
		double ta = (box0 - start) / dir;
		double tb = (box1 - start) / dir;
		if (ta > tb) {
			std::swap(ta, tb);
		}
		t0 = std::max(t0, ta);
		t1 = std::min(t1, tb);
		if (t0 > t1) {
			return false;
		}
	}
	return true;
}

// Distance to the piece axis, then the two end planes. A missing neighbor
// contributes a zero tangent, so that plane is perpendicular to the piece.
bool pointInEdelvikTube(
	const double* x,
	const double* a,
	const double* b,
	const double* tangent,
	double length,
	double radius,
	double tube_radius,
	const double* normal_left,
	const double* normal_right,
	double& rho,
	double& axial)
{
	double offset[3] = {x[0] - a[0], x[1] - a[1], x[2] - a[2]};
	const double s = dot3(offset, tangent);
	double radial2 = dot3(offset, offset) - s * s;
	if (radial2 < 0.0) {
		radial2 = 0.0;
	}
	rho = std::sqrt(radial2);
	if (rho < radius || rho > tube_radius) {
		return false;
	}
	const double left = (x[0] - a[0]) * normal_left[0]
		+ (x[1] - a[1]) * normal_left[1]
		+ (x[2] - a[2]) * normal_left[2];
	if (left < 0.0) {
		return false;
	}
	const double right = (x[0] - b[0]) * normal_right[0]
		+ (x[1] - b[1]) * normal_right[1]
		+ (x[2] - b[2]) * normal_right[2];
	if (right > 0.0) {
		return false;
	}
	axial = std::min(length, std::max(0.0, s));
	return true;
}

}  // namespace

double edelvikNormalizedInductance(double radius, double transverse)
{
	if (!(radius > 0.0) || !(transverse > radius)) {
		throw std::runtime_error(
			"Thin-wire radius must be positive and smaller than the transverse length.");
	}
	const double tube = 1.7 * transverse;
	if (!(tube > radius)) {
		throw std::runtime_error("Edelvik tube radius is not larger than the wire radius.");
	}
	const double inductance = std::log((tube + radius) / (2.0 * radius)) / (2.0 * M_PI);
	if (!(inductance > 0.0)) {
		throw std::runtime_error("Computed thin-wire inductance is not positive.");
	}
	return inductance;
}

double edelvikTubeWeight(double rho, double radius, double tube_radius)
{
	if (!(rho >= radius) || rho > tube_radius || !(tube_radius > radius)) {
		return 0.0;
	}
	const double theta = M_PI * radius / tube_radius;
	const double r2 = tube_radius * tube_radius;
	const double denom =
		M_PI * (r2 - radius * radius)
		- (2.0 * r2 / M_PI) * (1.0 + std::cos(theta) + theta * std::sin(theta));
	if (!(denom > 0.0)) {
		return 0.0;
	}
	return (1.0 + std::cos(M_PI * rho / tube_radius)) / denom;
}

std::vector<WirePolyline> buildWirePolylinesFromGmsh(
	const std::string& mesh_path,
	const std::vector<WireSpec>& specs)
{
	std::unordered_map<int, std::array<double, 3>> nodes;
	std::vector<CurveEdge> edges;
	readGmshCurves(mesh_path, nodes, edges);
	std::set<int> seen;
	std::vector<WirePolyline> curves;
	curves.reserve(specs.size());
	for (const auto& spec : specs) {
		for (int tag : spec.tags) {
			if (!seen.insert(tag).second) {
				throw std::runtime_error("A thin-wire elementary tag is listed twice.");
			}
		}
		curves.push_back(chainEdges(spec.tags, edges, nodes));
	}
	return curves;
}

std::vector<WirePieceGeom> buildWirePieceChain(const Model& model)
{
	std::vector<WirePieceGeom> pieces;
	if (model.wireSpecs().empty()) {
		return pieces;
	}
	const mfem::Mesh& mesh = model.getConstSerialMesh();
	if (mesh.Dimension() != 3) {
		throw std::runtime_error("Thin wires require a 3D volume mesh.");
	}
	const MPI_Comm comm = model.getConstMesh().GetComm();
	int comm_size = 1;
	MPI_Comm_size(comm, &comm_size);

	std::unordered_set<int> blocked;
	auto block_marker = [&](const mfem::Array<int>& marker) {
		const mfem::Array<int> united = uniteMarker(marker, comm, comm_size);
		if (united.Size() == 0) {
			return;
		}
		for (int be = 0; be < mesh.GetNBE(); ++be) {
			const int attr = mesh.GetBdrAttribute(be);
			if (!markerHits(united, attr)) {
				continue;
			}
			int elem1 = -1;
			int elem2 = -1;
			mesh.GetFaceElements(mesh.GetBdrElementFaceIndex(be), &elem1, &elem2);
			if (elem1 >= 0) {
				blocked.insert(elem1);
			}
			if (elem2 >= 0) {
				blocked.insert(elem2);
			}
		}
	};
	for (const auto& entry : model.getTotalFieldScatteredFieldToMarker()) {
		block_marker(entry.second);
	}
	for (const auto& entry : model.getSGBCToMarker()) {
		block_marker(entry.second);
	}
	block_marker(model.getDeltaGapMarker());

	struct RawPiece {
		double a[3];
		double b[3];
		int host = -1;
		double edge = 0.0;
		bool left_end = false;
		bool right_end = false;
		int wire = 0;
	};
	std::vector<RawPiece> raw;
	for (int w = 0; w < static_cast<int>(model.wireCurves().size()); ++w) {
		const auto& curve = model.wireCurves()[w];
		const int nseg = static_cast<int>(curve.points.size()) - 1;
		for (int s = 0; s < nseg; ++s) {
			const double* a = curve.points[s].data();
			const double* b = curve.points[s + 1].data();
			std::vector<Clip> hits;
			for (int el = 0; el < mesh.GetNE(); ++el) {
				mfem::Array<int> verts;
				mesh.GetElementVertices(el, verts);
				if (verts.Size() != 4) {
					throw std::runtime_error("Thin wires require a tetrahedral mesh.");
				}
				double coordinates[4][3];
				for (int v = 0; v < 4; ++v) {
					const double* coord = mesh.GetVertex(verts[v]);
					coordinates[v][0] = coord[0];
					coordinates[v][1] = coord[1];
					coordinates[v][2] = coord[2];
				}
				double s0 = 0.0;
				double s1 = 1.0;
				if (!clipSegmentToTet(a, b, coordinates, s0, s1)) {
					continue;
				}
				if (s1 - s0 < 1e-10) {
					continue;
				}
				hits.push_back(Clip{el, std::max(0.0, s0), std::min(1.0, s1)});
			}
			if (hits.empty()) {
				throw std::runtime_error("A thin-wire segment does not cross any tetrahedron.");
			}
			std::sort(hits.begin(), hits.end(), [](const Clip& lhs, const Clip& rhs) {
				return lhs.s0 < rhs.s0;
			});
			double cursor = 0.0;
			for (const Clip& hit : hits) {
				const double start = std::max(cursor, hit.s0);
				const double stop = hit.s1;
				if (stop - start < 1e-10) {
					continue;
				}
				if (blocked.count(hit.element)) {
					throw std::runtime_error(
						"A thin wire crosses a TFSF, SGBC, or delta-gap element.");
				}
				const int attr = mesh.GetAttribute(hit.element);
				if (model.findDebye(attr) || model.findLorentz(attr)) {
					throw std::runtime_error("A thin wire must stay outside Debye and Lorentz materials.");
				}
				if (const PMLProperties* pml = model.getPMLPropertiesForTag(attr)) {
					if (std::abs(pml->kappa_max - 1.0) > 1e-12) {
						throw std::runtime_error(
							"A thin wire in PML requires kappa_max = 1. "
							"Above 1 the electric update is divided by the stretch, and the wire current is not.");
					}
				} else {
					const auto& materials = model.getGeomTagToMaterial();
					const auto material = materials.find(attr);
					if (material != materials.end()) {
						const Material& med = material->second;
						if (med.getPermittivity() != 1.0 || med.getPermeability() != 1.0
							|| med.getConductivity() != 0.0) {
							throw std::runtime_error("A thin wire must stay in vacuum.");
						}
					}
				}
				RawPiece piece;
				for (int d = 0; d < 3; ++d) {
					piece.a[d] = a[d] + start * (b[d] - a[d]);
					piece.b[d] = a[d] + stop * (b[d] - a[d]);
				}
				piece.host = hit.element;
				mfem::Array<int> host_verts;
				mesh.GetElementVertices(hit.element, host_verts);
				double coordinates[4][3];
				for (int v = 0; v < 4; ++v) {
					const double* coord = mesh.GetVertex(host_verts[v]);
					coordinates[v][0] = coord[0];
					coordinates[v][1] = coord[1];
					coordinates[v][2] = coord[2];
				}
				piece.edge = meanEdgeLength(coordinates);
				piece.wire = w;
				raw.push_back(piece);
				cursor = stop;
			}
			if (cursor < 0.99) {
				throw std::runtime_error("A thin-wire segment is not fully covered by the mesh.");
			}
		}
	}
	if (raw.empty()) {
		throw std::runtime_error("The thin-wire curve produced no segments inside the mesh.");
	}

	std::vector<int> wire_begin(model.wireCurves().size(), -1);
	std::vector<int> wire_end(model.wireCurves().size(), -1);
	for (int i = 0; i < static_cast<int>(raw.size()); ++i) {
		if (wire_begin[raw[i].wire] < 0) {
			wire_begin[raw[i].wire] = i;
		}
		wire_end[raw[i].wire] = i;
	}
	for (int w = 0; w < static_cast<int>(model.wireCurves().size()); ++w) {
		if (wire_begin[w] < 0) {
			continue;
		}
		if (!model.wireCurves()[w].closed) {
			raw[wire_begin[w]].left_end = true;
			raw[wire_end[w]].right_end = true;
		}
	}

	pieces.resize(raw.size());
	for (int p = 0; p < static_cast<int>(raw.size()); ++p) {
		WirePieceGeom& piece = pieces[p];
		const RawPiece& source = raw[p];
		piece.host = source.host;
		piece.wire = source.wire;
		const int begin = wire_begin[source.wire];
		const int end = wire_end[source.wire];
		const bool closed = model.wireCurves()[source.wire].closed;
		if (source.left_end) {
			piece.left_neighbor = -1;
		} else if (closed && p == begin) {
			piece.left_neighbor = end;
		} else {
			piece.left_neighbor = p - 1;
		}
		if (source.right_end) {
			piece.right_neighbor = -1;
		} else if (closed && p == end) {
			piece.right_neighbor = begin;
		} else {
			piece.right_neighbor = p + 1;
		}
		double delta[3] = {
			source.b[0] - source.a[0],
			source.b[1] - source.a[1],
			source.b[2] - source.a[2]};
		piece.length = norm3(delta);
		if (!(piece.length > 0.0)) {
			throw std::runtime_error("A thin-wire piece has zero length.");
		}
		for (int d = 0; d < 3; ++d) {
			piece.tangent[d] = delta[d] / piece.length;
			piece.a[d] = source.a[d];
			piece.b[d] = source.b[d];
		}
		const WireSpec& spec = model.wireSpecs()[source.wire];
		const double cell = spec.transverse_step > 0.0 ? spec.transverse_step : source.edge;
		const double geometric = edelvikNormalizedInductance(spec.radius, cell);
		const double extra = spec.inductance_per_meter / physicalConstants::vacuumPermeability_SI;
		piece.resistance = spec.resistance_per_meter / physicalConstants::freeSpaceImpedance_SI;
		piece.radius = spec.radius;
		piece.tube_radius = 1.7 * cell;
		const double total = geometric + extra;
		piece.v2 = geometric / total;
		piece.inductance = total;
		if (source.left_end) {
			piece.left_pec = endCondition(mesh, model, source.a) == 1;
		}
		if (source.right_end) {
			piece.right_pec = endCondition(mesh, model, source.b) == 1;
		}
	}
	return pieces;
}

ThinWireCoupler::ThinWireCoupler(mfem::ParFiniteElementSpace& fes, const Model& model)
	: model_(model)
{
	if (model.wireSpecs().empty()) {
		return;
	}
	if (fes.GetMesh()->Dimension() != 3) {
		throw std::runtime_error("Thin wires require a 3D volume mesh.");
	}
	const auto* collection = dynamic_cast<const mfem::DG_FECollection*>(fes.FEColl());
	if (collection == nullptr) {
		throw std::runtime_error("Thin wires require a DG finite-element space.");
	}
	order_ = collection->GetOrder();
	basis_type_ = collection->GetBasisType();
	if (order_ < 1) {
		throw std::runtime_error("Thin wires require polynomial order at least 1.");
	}
	comm_ = fes.GetParMesh()->GetComm();
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (mfem::Device::Allows(mfem::Backend::CUDA)) {
		int comm_size = 1;
		MPI_Comm_size(comm_, &comm_size);
		if (comm_size > 1) {
			throw std::runtime_error("Thin wires on CUDA require a single MPI rank.");
		}
		device_path_ = true;
	}
#endif

	std::vector<WirePieceGeom> local_chain;
	const std::vector<WirePieceGeom>* chain = &model.wirePieces();
	if (chain->empty()) {
		local_chain = buildWirePieceChain(model);
		chain = &local_chain;
	}

	mfem::DG_FECollection line_collection(order_, 1, collection->GetBasisType());
	mfem::ParMesh& mesh = *fes.GetParMesh();
	mfem::Mesh line_mesh(1, static_cast<int>(chain->size()) + 1, static_cast<int>(chain->size()), 0, 3);
	for (int i = 0; i < static_cast<int>(chain->size()); ++i) {
		if (i == 0) {
			line_mesh.AddVertex((*chain)[i].a);
		}
		line_mesh.AddVertex((*chain)[i].b);
		line_mesh.AddSegment(i, i + 1, 1);
	}
	line_mesh.FinalizeTopology();
	mfem::FiniteElementSpace line_space(&line_mesh, &line_collection);
	const int nloc = order_ + 1;
	const mfem::FiniteElement& line_fe = *line_space.GetFE(0);
	mfem::IntegrationPoint left_ip;
	mfem::IntegrationPoint right_ip;
	left_ip.Set1w(0.0, 0.0);
	right_ip.Set1w(1.0, 0.0);
	mfem::Vector shape_left(nloc);
	mfem::Vector shape_right(nloc);
	line_fe.CalcShape(left_ip, shape_left);
	line_fe.CalcShape(right_ip, shape_right);
	const mfem::IntegrationRule& line_rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * order_ + 1);
	const mfem::IntegrationRule& vol_rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, 2 * order_);
	const mfem::IntegrationRule& rho_rule = mfem::IntRules.Get(
		mfem::Geometry::SEGMENT, std::max(2 * order_ + 2, 12));
	const int n_theta = std::max(16, 8 * (order_ + 1));
	const double theta_step = 2.0 * std::acos(-1.0) / static_cast<double>(n_theta);

	pieces_.resize(chain->size());
	int offset = 0;
	double min_length = 1e300;
	for (int p = 0; p < static_cast<int>(chain->size()); ++p) {
		Piece& piece = pieces_[p];
		const WirePieceGeom& source = (*chain)[p];
		const int local_host = serialElementToLocal(model, mesh, source.host);
		piece.owned = local_host >= 0;
		piece.host = local_host;
		piece.wire = source.wire;
		piece.left_neighbor = source.left_neighbor;
		piece.right_neighbor = source.right_neighbor;
		piece.left_pec = source.left_pec;
		piece.right_pec = source.right_pec;
		piece.length = source.length;
		piece.inductance = source.inductance;
		piece.resistance = source.resistance;
		piece.v2 = source.v2;
		piece.radius = source.radius;
		piece.tube_radius = source.tube_radius;
		min_length = std::min(min_length, piece.length);
		for (int d = 0; d < 3; ++d) {
			piece.tangent[d] = source.tangent[d];
			piece.a[d] = source.a[d];
			piece.b[d] = source.b[d];
		}
		piece.offset = offset;
		offset += nloc;
		piece.shape_left = shape_left;
		piece.shape_right = shape_right;

		mfem::DenseMatrix line_mass(nloc);
		mfem::DenseMatrix line_stiff(nloc);
		line_mass = 0.0;
		line_stiff = 0.0;
		mfem::ElementTransformation* line_trans = line_mesh.GetElementTransformation(p);
		for (int q = 0; q < line_rule.GetNPoints(); ++q) {
			const mfem::IntegrationPoint& ip = line_rule.IntPoint(q);
			mfem::Vector shape(nloc);
			mfem::DenseMatrix dshape(nloc, 1);
			line_fe.CalcShape(ip, shape);
			line_fe.CalcDShape(ip, dshape);
			line_trans->SetIntPoint(&ip);
			const double weight = ip.weight;
			for (int i = 0; i < nloc; ++i) {
				for (int j = 0; j < nloc; ++j) {
					line_mass(i, j) += shape(i) * shape(j) * piece.length * weight;
					line_stiff(i, j) += shape(i) * dshape(j, 0) * weight;
				}
			}
		}
		piece.mass_inv.SetSize(nloc);
		piece.mass_inv = line_mass;
		piece.mass_inv.Invert();
		piece.stiffness = line_stiff;
	}

	const int nelem = mesh.GetNE();
	std::vector<mfem::DenseMatrix> vol_mass_inv(static_cast<std::size_t>(nelem));
	std::vector<char> have_mass(static_cast<std::size_t>(nelem), 0);
	int tube_blocks = 0;
	for (int p = 0; p < static_cast<int>(pieces_.size()); ++p) {
		Piece& piece = pieces_[static_cast<std::size_t>(p)];
		double normal_left[3] = {piece.tangent[0], piece.tangent[1], piece.tangent[2]};
		double normal_right[3] = {piece.tangent[0], piece.tangent[1], piece.tangent[2]};
		if (piece.left_neighbor >= 0) {
			const Piece& other = pieces_[static_cast<std::size_t>(piece.left_neighbor)];
			for (int d = 0; d < 3; ++d) {
				normal_left[d] += other.tangent[d];
			}
		}
		if (piece.right_neighbor >= 0) {
			const Piece& other = pieces_[static_cast<std::size_t>(piece.right_neighbor)];
			for (int d = 0; d < 3; ++d) {
				normal_right[d] += other.tangent[d];
			}
		}
		double axis_hint[3] = {0.0, 0.0, 1.0};
		if (std::abs(piece.tangent[2]) > 0.9) {
			axis_hint[2] = 0.0;
			axis_hint[0] = 1.0;
		}
		double normal_seed[3];
		cross3(piece.tangent, axis_hint, normal_seed);
		const double seed_norm = norm3(normal_seed);
		if (!(seed_norm > 0.0)) {
			throw std::runtime_error("Thin-wire tube frame is undefined.");
		}
		double radial_u[3];
		double radial_v[3];
		for (int d = 0; d < 3; ++d) {
			radial_u[d] = normal_seed[d] / seed_norm;
		}
		cross3(piece.tangent, radial_u, radial_v);

		std::vector<int> candidates;
		for (int elem = 0; elem < nelem; ++elem) {
			mfem::Array<int> vertices;
			mesh.GetElementVertices(elem, vertices);
			double box_lo[3] = {1e300, 1e300, 1e300};
			double box_hi[3] = {-1e300, -1e300, -1e300};
			for (int v = 0; v < vertices.Size(); ++v) {
				const double* coord = mesh.GetVertex(vertices[v]);
				for (int d = 0; d < 3; ++d) {
					box_lo[d] = std::min(box_lo[d], coord[d]);
					box_hi[d] = std::max(box_hi[d], coord[d]);
				}
			}
			if (segmentMeetsExpandedBox(piece.a, piece.b, box_lo, box_hi, piece.tube_radius)) {
				candidates.push_back(elem);
			}
		}
		if (candidates.empty()) {
			continue;
		}

		struct TubeAccum {
			int elem = -1;
			mfem::Array<int> dofs;
			mfem::DenseMatrix load_x, load_y, load_z;
			mfem::DenseMatrix sample_x, sample_y, sample_z;
			bool hit = false;
		};
		std::vector<TubeAccum> accum(candidates.size());
		for (std::size_t c = 0; c < candidates.size(); ++c) {
			const int elem = candidates[c];
			TubeAccum& slot = accum[c];
			slot.elem = elem;
			fes.GetElementDofs(elem, slot.dofs);
			for (int i = 0; i < slot.dofs.Size(); ++i) {
				const int dof = slot.dofs[i];
				slot.dofs[i] = dof >= 0 ? dof : -1 - dof;
			}
			const int nvol = slot.dofs.Size();
			slot.load_x.SetSize(nvol, nloc);
			slot.load_y.SetSize(nvol, nloc);
			slot.load_z.SetSize(nvol, nloc);
			slot.sample_x.SetSize(nloc, nvol);
			slot.sample_y.SetSize(nloc, nvol);
			slot.sample_z.SetSize(nloc, nvol);
			slot.load_x = 0.0;
			slot.load_y = 0.0;
			slot.load_z = 0.0;
			slot.sample_x = 0.0;
			slot.sample_y = 0.0;
			slot.sample_z = 0.0;
		}

		const double radial_span = piece.tube_radius - piece.radius;
		for (int qs = 0; qs < line_rule.GetNPoints(); ++qs) {
			const mfem::IntegrationPoint& sip = line_rule.IntPoint(qs);
			const double axial = sip.x * piece.length;
			const double ds = sip.weight * piece.length;
			mfem::Vector lshape(nloc);
			line_fe.CalcShape(sip, lshape);
			for (int it = 0; it < n_theta; ++it) {
				const double theta = theta_step * static_cast<double>(it);
				const double ct = std::cos(theta);
				const double st = std::sin(theta);
				for (int qr = 0; qr < rho_rule.GetNPoints(); ++qr) {
					const mfem::IntegrationPoint& rip = rho_rule.IntPoint(qr);
					const double rho = piece.radius + rip.x * radial_span;
					const double dr = rip.weight * radial_span;
					const double gain = edelvikTubeWeight(rho, piece.radius, piece.tube_radius);
					if (!(gain > 0.0)) {
						continue;
					}
					double phys[3];
					for (int d = 0; d < 3; ++d) {
						phys[d] = piece.a[d] + piece.tangent[d] * axial
							+ (ct * radial_u[d] + st * radial_v[d]) * rho;
					}
					double rho_check = 0.0;
					double axial_check = 0.0;
					if (!pointInEdelvikTube(
							phys, piece.a, piece.b, piece.tangent, piece.length,
							piece.radius, piece.tube_radius, normal_left, normal_right,
							rho_check, axial_check)) {
						continue;
					}
					const double measure = gain * rho * dr * theta_step * ds;
					int found = -1;
					mfem::IntegrationPoint vol_ip;
					for (std::size_t c = 0; c < candidates.size(); ++c) {
						auto* trans = mesh.GetElementTransformation(candidates[c]);
						mfem::InverseElementTransformation inverse;
						inverse.SetTransformation(*trans);
						mfem::Vector point(3);
						point(0) = phys[0];
						point(1) = phys[1];
						point(2) = phys[2];
						if (inverse.Transform(point, vol_ip)
							== mfem::InverseElementTransformation::Inside) {
							found = static_cast<int>(c);
							break;
						}
					}
					if (found < 0) {
						continue;
					}
					TubeAccum& slot = accum[static_cast<std::size_t>(found)];
					const int nvol = slot.dofs.Size();
					const mfem::FiniteElement& vol_fe = *fes.GetFE(slot.elem);
					mfem::Vector vshape(nvol);
					vol_fe.CalcShape(vol_ip, vshape);
					slot.hit = true;
					for (int i = 0; i < nvol; ++i) {
						for (int j = 0; j < nloc; ++j) {
							const double value = vshape(i) * lshape(j) * measure;
							slot.load_x(i, j) += piece.tangent[0] * value;
							slot.sample_x(j, i) += piece.tangent[0] * value;
							slot.load_y(i, j) += piece.tangent[1] * value;
							slot.sample_y(j, i) += piece.tangent[1] * value;
							slot.load_z(i, j) += piece.tangent[2] * value;
							slot.sample_z(j, i) += piece.tangent[2] * value;
						}
					}
				}
			}
		}

		for (TubeAccum& slot : accum) {
			if (!slot.hit) {
				continue;
			}
			const int nvol = slot.dofs.Size();
			if (!have_mass[static_cast<std::size_t>(slot.elem)]) {
				mfem::DenseMatrix mass(nvol);
				mass = 0.0;
				const mfem::FiniteElement& vol_fe = *fes.GetFE(slot.elem);
				auto* mass_trans = mesh.GetElementTransformation(slot.elem);
				for (int q = 0; q < vol_rule.GetNPoints(); ++q) {
					const mfem::IntegrationPoint& ip = vol_rule.IntPoint(q);
					mass_trans->SetIntPoint(&ip);
					mfem::Vector shape(nvol);
					vol_fe.CalcShape(ip, shape);
					const double weight = ip.weight * mass_trans->Weight();
					for (int i = 0; i < nvol; ++i) {
						for (int j = 0; j < nvol; ++j) {
							mass(i, j) += shape(i) * shape(j) * weight;
						}
					}
				}
				mass.Invert();
				vol_mass_inv[static_cast<std::size_t>(slot.elem)] = mass;
				have_mass[static_cast<std::size_t>(slot.elem)] = 1;
			}
			TubeBlock block;
			block.vol_dofs = slot.dofs;
			block.es_x = slot.sample_x;
			block.es_y = slot.sample_y;
			block.es_z = slot.sample_z;
			mfem::DenseMatrix mixed(nvol, nloc);
			const mfem::DenseMatrix& mass_inv_vol = vol_mass_inv[static_cast<std::size_t>(slot.elem)];
			mfem::Mult(mass_inv_vol, slot.load_x, mixed);
			block.j_x = mixed;
			block.j_x *= -1.0;
			mfem::Mult(mass_inv_vol, slot.load_y, mixed);
			block.j_y = mixed;
			block.j_y *= -1.0;
			mfem::Mult(mass_inv_vol, slot.load_z, mixed);
			block.j_z = mixed;
			block.j_z *= -1.0;
			piece.blocks.push_back(std::move(block));
			++tube_blocks;
		}
	}

	if (mfem::Mpi::WorldRank() == 0) {
		std::cout << "Thin wires: " << pieces_.size()
		          << " pieces, " << offset << " line degrees of freedom, "
		          << tube_blocks << " tube elements." << std::endl;
	}

	current_.SetSize(offset);
	charge_.SetSize(offset);
	stage_current_.SetSize(offset);
	stage_charge_.SetSize(offset);
	d_current_.SetSize(offset);
	d_charge_.SetSize(offset);
	if (device_path_) {
		current_.UseDevice(true);
		charge_.UseDevice(true);
		stage_current_.UseDevice(true);
		stage_charge_.UseDevice(true);
		d_current_.UseDevice(true);
		d_charge_.UseDevice(true);
	}
	current_ = 0.0;
	charge_ = 0.0;
	const double speed = 1.0;
	substep_ = 0.5 * min_length / (speed * static_cast<double>((order_ + 1) * (order_ + 1)));
	if (device_path_) {
		packDevice();
	}
	buildExportMeshes();
}

void ThinWireCoupler::testingSetElectricReturnScale(double scale)
{
	electric_return_scale_ = scale;
}

void ThinWireCoupler::testingSetUniformCurrent(double value)
{
	double* current = current_.HostWrite();
	double* charge = charge_.HostWrite();
	for (int i = 0; i < current_.Size(); ++i) {
		current[i] = value;
		charge[i] = 0.0;
	}
}

void ThinWireCoupler::testingSetOpenHalfWave()
{
	double total = 0.0;
	for (const Piece& piece : pieces_) {
		total += piece.length;
	}
	if (!(total > 0.0)) {
		throw std::runtime_error("A half-wave needs a wire of positive length.");
	}
	const int nloc = order_ + 1;
	mfem::DG_FECollection collection(order_, 1, basis_type_);
	mfem::Mesh segment(1, 2, 1);
	segment.AddVertex(0.0);
	segment.AddVertex(1.0);
	segment.AddSegment(0, 1, 1);
	segment.FinalizeTopology();
	mfem::FiniteElementSpace space(&segment, &collection);
	const mfem::IntegrationRule& nodes = space.GetFE(0)->GetNodes();
	double* current = current_.HostWrite();
	double* charge = charge_.HostWrite();
	for (int i = 0; i < charge_.Size(); ++i) {
		charge[i] = 0.0;
	}
	double distance = 0.0;
	for (const Piece& piece : pieces_) {
		for (int i = 0; i < nloc; ++i) {
			const double xi = nodes.IntPoint(i).x;
			const double s = distance + xi * piece.length;
			current[piece.offset + i] = std::sin(M_PI * s / total);
		}
		distance += piece.length;
	}
}

double ThinWireCoupler::testingPieceLength(int piece) const
{
	return pieces_.at(static_cast<std::size_t>(piece)).length;
}

double ThinWireCoupler::testingPieceInductance(int piece) const
{
	return pieces_.at(static_cast<std::size_t>(piece)).inductance;
}

int ThinWireCoupler::testingOwnedCount() const
{
	int count = 0;
	for (const Piece& piece : pieces_) {
		if (piece.owned) {
			++count;
		}
	}
	return count;
}

double ThinWireCoupler::testingCurrentValue(int dof) const
{
	return current_.HostRead()[dof];
}

double ThinWireCoupler::testingStageCurrentValue(int dof) const
{
	return stage_current_.HostRead()[dof];
}

double ThinWireCoupler::testingStageChargeValue(int dof) const
{
	return stage_charge_.HostRead()[dof];
}

double ThinWireCoupler::testingCurrentNorm() const
{
	return current_.Norml2();
}

void ThinWireCoupler::addCurrentToElectric(
	const mfem::Vector& current, mfem::Vector& out, int ndofs) const
{
	for (const Piece& piece : pieces_) {
		for (const TubeBlock& block : piece.blocks) {
			const mfem::DenseMatrix* matrices[3] = {&block.j_x, &block.j_y, &block.j_z};
			for (int d = 0; d < 3; ++d) {
				for (int i = 0; i < block.vol_dofs.Size(); ++i) {
					double sum = 0.0;
					for (int j = 0; j < matrices[d]->Width(); ++j) {
						sum += (*matrices[d])(i, j) * current(piece.offset + j);
					}
					const int dof = block.vol_dofs[i];
					out[d * ndofs + dof] += electric_return_scale_ * sum;
				}
			}
		}
	}
}

void ThinWireCoupler::electricDrive(
	const mfem::Vector* eh, int ndofs, mfem::Vector& drive) const
{
	drive.SetSize(current_.Size());
	drive = 0.0;
	if (eh == nullptr || current_.Size() == 0) {
		return;
	}
	const int nloc = order_ + 1;
	for (const Piece& piece : pieces_) {
		mfem::Vector e_load(nloc);
		e_load = 0.0;
		for (const TubeBlock& block : piece.blocks) {
			const mfem::DenseMatrix* sample[3] = {&block.es_x, &block.es_y, &block.es_z};
			for (int d = 0; d < 3; ++d) {
				for (int j = 0; j < nloc; ++j) {
					double sum = 0.0;
					for (int i = 0; i < block.vol_dofs.Size(); ++i) {
						sum += (*sample[d])(j, i) * (*eh)[d * ndofs + block.vol_dofs[i]];
					}
					e_load(j) += sum;
				}
			}
		}
		mfem::Vector e_nodal(nloc);
		piece.mass_inv.Mult(e_load, e_nodal);
		for (int i = 0; i < nloc; ++i) {
			drive(piece.offset + i) = e_nodal(i);
		}
	}
	int comm_size = 1;
	MPI_Comm_size(comm_, &comm_size);
	if (comm_size > 1) {
		MPI_Allreduce(MPI_IN_PLACE, drive.GetData(), drive.Size(), MPI_DOUBLE, MPI_SUM, comm_);
	}
}

void ThinWireCoupler::rhs(
	const mfem::Vector& current,
	const mfem::Vector& charge,
	const mfem::Vector* drive,
	mfem::Vector& d_current,
	mfem::Vector& d_charge) const
{
	d_current = 0.0;
	d_charge = 0.0;
	const int nloc = order_ + 1;
	auto local = [&](const mfem::Vector& src, const Piece& piece, mfem::Vector& dst) {
		dst.SetSize(nloc);
		for (int i = 0; i < nloc; ++i) {
			dst(i) = src(piece.offset + i);
		}
	};
	auto trace = [](const mfem::Vector& dofs, const mfem::Vector& shape) {
		double value = 0.0;
		for (int i = 0; i < dofs.Size(); ++i) {
			value += dofs(i) * shape(i);
		}
		return value;
	};

	for (int p = 0; p < static_cast<int>(pieces_.size()); ++p) {
		const Piece& piece = pieces_[p];
		mfem::Vector I(nloc);
		mfem::Vector q(nloc);
		local(current, piece, I);
		local(charge, piece, q);
		const double speed = std::sqrt(piece.v2);
		mfem::Vector stiff_q(nloc);
		mfem::Vector stiff_I(nloc);
		piece.stiffness.Mult(q, stiff_q);
		piece.stiffness.Mult(I, stiff_I);
		mfem::Vector vol_q(nloc);
		mfem::Vector vol_I(nloc);
		piece.mass_inv.Mult(stiff_q, vol_q);
		piece.mass_inv.Mult(stiff_I, vol_I);

		mfem::Vector dI(nloc);
		mfem::Vector dq(nloc);
		for (int i = 0; i < nloc; ++i) {
			dI(i) = -piece.v2 * vol_q(i) - (piece.resistance / piece.inductance) * I(i);
			dq(i) = -vol_I(i);
		}

		if (drive != nullptr) {
			for (int i = 0; i < nloc; ++i) {
				dI(i) += (*drive)(piece.offset + i) / piece.inductance;
			}
		}

		auto hat = [&](bool right, double& I_hat, double& q_hat) {
			const double I_i = trace(I, right ? piece.shape_right : piece.shape_left);
			const double q_i = trace(q, right ? piece.shape_right : piece.shape_left);
			const int neighbor = right ? piece.right_neighbor : piece.left_neighbor;
			const bool pec = right ? piece.right_pec : piece.left_pec;
			const double alpha_out = right
				? 0.5 * (q_i + I_i / speed)
				: 0.5 * (q_i - I_i / speed);
			double alpha_in = 0.0;
			if (neighbor >= 0) {
				const Piece& other = pieces_[neighbor];
				mfem::Vector I_n(nloc);
				mfem::Vector q_n(nloc);
				local(current, other, I_n);
				local(charge, other, q_n);
				const double I_e = trace(I_n, right ? other.shape_left : other.shape_right);
				const double q_e = trace(q_n, right ? other.shape_left : other.shape_right);
				alpha_in = right
					? 0.5 * (q_e - I_e / speed)
					: 0.5 * (q_e + I_e / speed);
			} else if (pec) {
				alpha_in = -alpha_out;
			} else {
				alpha_in = alpha_out;
			}
			if (right) {
				I_hat = speed * (alpha_out - alpha_in);
				q_hat = alpha_out + alpha_in;
			} else {
				I_hat = speed * (alpha_in - alpha_out);
				q_hat = alpha_in + alpha_out;
			}
		};

		auto add_end = [&](bool right) {
			double I_hat = 0.0;
			double q_hat = 0.0;
			hat(right, I_hat, q_hat);
			const mfem::Vector& shape = right ? piece.shape_right : piece.shape_left;
			const double I_i = trace(I, shape);
			const double q_i = trace(q, shape);
			const double normal = right ? 1.0 : -1.0;
			const double flux_I = normal * piece.v2 * (q_i - q_hat);
			const double flux_q = normal * (I_i - I_hat);
			mfem::Vector lift_I(nloc);
			mfem::Vector lift_q(nloc);
			for (int i = 0; i < nloc; ++i) {
				lift_I(i) = shape(i) * flux_I;
				lift_q(i) = shape(i) * flux_q;
			}
			mfem::Vector add_I(nloc);
			mfem::Vector add_q(nloc);
			piece.mass_inv.Mult(lift_I, add_I);
			piece.mass_inv.Mult(lift_q, add_q);
			dI.Add(1.0, add_I);
			dq.Add(1.0, add_q);
		};
		add_end(false);
		add_end(true);

		for (int i = 0; i < nloc; ++i) {
			d_current(piece.offset + i) = dI(i);
			d_charge(piece.offset + i) = dq(i);
		}
	}
}

double ThinWireCoupler::feedElectric(double time) const
{
	if (feed_count_ <= 0) {
		return 0.0;
	}
	return gaussianTimeSignal(
		time, feed_t0_, feed_spread_, feed_magnitude_, feed_derivative_);
}

void ThinWireCoupler::setFeed(
	const double position[3], double magnitude, double spread, double t0, bool derivative)
{
	if (!(magnitude > 0.0) || !std::isfinite(magnitude)) {
		throw std::runtime_error("wire_voltage magnitude must be > 0.");
	}
	if (!(spread > 0.0) || !std::isfinite(spread)) {
		throw std::runtime_error("wire_voltage spread must be > 0.");
	}
	if (!std::isfinite(t0)) {
		throw std::runtime_error("wire_voltage t0 must be finite.");
	}
	int hits[kMaxFeeds];
	int hit_count = 0;
	for (int p = 0; p < static_cast<int>(pieces_.size()); ++p) {
		const Piece& piece = pieces_[p];
		double ab[3] = {
			piece.b[0] - piece.a[0],
			piece.b[1] - piece.a[1],
			piece.b[2] - piece.a[2]};
		const double length2 = dot3(ab, ab);
		if (!(length2 > 0.0)) {
			continue;
		}
		double ap[3] = {
			position[0] - piece.a[0],
			position[1] - piece.a[1],
			position[2] - piece.a[2]};
		const double xi = dot3(ap, ab) / length2;
		if (xi < -1e-8 || xi > 1.0 + 1e-8) {
			continue;
		}
		const double clamped = std::max(0.0, std::min(1.0, xi));
		double closest[3] = {
			piece.a[0] + clamped * ab[0],
			piece.a[1] + clamped * ab[1],
			piece.a[2] + clamped * ab[2]};
		double delta[3] = {
			position[0] - closest[0],
			position[1] - closest[1],
			position[2] - closest[2]};
		if (norm3(delta) <= 1e-8) {
			if (hit_count >= kMaxFeeds) {
				throw std::runtime_error("wire_voltage position meets more than two pieces.");
			}
			hits[hit_count++] = p;
		}
	}
	if (hit_count == 0) {
		throw std::runtime_error("wire_voltage position is not on a wire.");
	}
	feed_count_ = hit_count;
	const double share = 1.0 / static_cast<double>(hit_count);
	for (int i = 0; i < hit_count; ++i) {
		feed_piece_[i] = hits[i];
		feed_per_length_[i] = share / pieces_[hits[i]].length;
	}
	for (int i = hit_count; i < kMaxFeeds; ++i) {
		feed_piece_[i] = -1;
		feed_per_length_[i] = 0.0;
	}
	feed_magnitude_ = magnitude;
	feed_spread_ = spread;
	feed_t0_ = t0;
	feed_derivative_ = derivative;
	if (mfem::Mpi::WorldRank() == 0) {
		std::cout << "[wire_voltage] pieces=" << hit_count
		          << " first=" << feed_piece_[0]
		          << " length=" << pieces_[feed_piece_[0]].length
		          << " magnitude=" << magnitude
		          << " spread=" << spread
		          << " t0=" << t0 << "\n";
	}
}

void ThinWireCoupler::advance(
	double tau,
	const mfem::Vector* eh,
	int ndofs,
	mfem::Vector& current,
	mfem::Vector& charge,
	double sample_time) const
{
	if (!(tau > 0.0)) {
		return;
	}
	int steps = std::max(1, static_cast<int>(std::ceil(tau / substep_)));
	const double h = tau / static_cast<double>(steps);
	mfem::Vector drive;
	electricDrive(eh, ndofs, drive);
	if (feed_count_ > 0) {
		const double voltage = feedElectric(sample_time);
		const int nloc = order_ + 1;
		for (int feed = 0; feed < feed_count_; ++feed) {
			const double extra = voltage * feed_per_length_[feed];
			const int off = pieces_[feed_piece_[feed]].offset;
			for (int i = 0; i < nloc; ++i) {
				drive(off + i) += extra;
			}
		}
	}
	for (int step = 0; step < steps; ++step) {
		rhs(current, charge, &drive, d_current_, d_charge_);
		current.Add(h, d_current_);
		charge.Add(h, d_charge_);
	}
}

void ThinWireCoupler::commit(double time, double dt)
{
	step_time_ = time;
	step_dt_ = dt;
}

void ThinWireCoupler::contributeHost(
	const mfem::Vector& in, mfem::Vector& out, int ndofs, double stage_time) const
{
	const double tau = std::max(0.0, stage_time - step_time_);
	in.HostRead();
	out.HostReadWrite();
	const bool current_dev = current_.UseDevice();
	const bool charge_dev = charge_.UseDevice();
	const bool stage_current_dev = stage_current_.UseDevice();
	const bool stage_charge_dev = stage_charge_.UseDevice();
	const bool d_current_dev = d_current_.UseDevice();
	const bool d_charge_dev = d_charge_.UseDevice();
	current_.HostRead();
	charge_.HostRead();
	current_.UseDevice(false);
	charge_.UseDevice(false);
	stage_current_.UseDevice(false);
	stage_charge_.UseDevice(false);
	d_current_.UseDevice(false);
	d_charge_.UseDevice(false);
	stage_current_ = current_;
	stage_charge_ = charge_;
	advance(tau, &in, ndofs, stage_current_, stage_charge_, stage_time);
	addCurrentToElectric(stage_current_, out, ndofs);
	current_.UseDevice(current_dev);
	charge_.UseDevice(charge_dev);
	stage_current_.UseDevice(stage_current_dev);
	stage_charge_.UseDevice(stage_charge_dev);
	d_current_.UseDevice(d_current_dev);
	d_charge_.UseDevice(d_charge_dev);
}

void ThinWireCoupler::contribute(
	const mfem::Vector& in, mfem::Vector& out, int ndofs, double stage_time) const
{
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (device_path_) {
		const double tau = std::max(0.0, stage_time - step_time_);
		advanceOnDevice(tau, in, &out, ndofs, stage_current_, stage_charge_, true, stage_time);
		return;
	}
#endif
	contributeHost(in, out, ndofs, stage_time);
}

void ThinWireCoupler::testingContributeOnHost(
	const mfem::Vector& in, mfem::Vector& out, int ndofs, double stage_time) const
{
	contributeHost(in, out, ndofs, stage_time);
}

void ThinWireCoupler::finalize(const mfem::Vector& eh, int ndofs)
{
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (device_path_) {
		advanceOnDevice(step_dt_, eh, nullptr, ndofs, current_, charge_, false,
		step_time_ + step_dt_);
		return;
	}
#endif
	advance(step_dt_, &eh, ndofs, current_, charge_, step_time_ + step_dt_);
}

bool ThinWireCoupler::sample(const std::vector<double>& point, double& current, double& charge) const
{
	if (point.size() < 3 || pieces_.empty()) {
		return false;
	}
	const double p[3] = {point[0], point[1], point[2]};
	int best = -1;
	double best_xi = 0.0;
	double best_distance = 1e300;
	for (int n = 0; n < static_cast<int>(pieces_.size()); ++n) {
		const Piece& piece = pieces_[n];
		double ab[3] = {
			piece.b[0] - piece.a[0],
			piece.b[1] - piece.a[1],
			piece.b[2] - piece.a[2]};
		const double length2 = dot3(ab, ab);
		if (!(length2 > 0.0)) {
			continue;
		}
		double ap[3] = {p[0] - piece.a[0], p[1] - piece.a[1], p[2] - piece.a[2]};
		double xi = dot3(ap, ab) / length2;
		xi = std::max(0.0, std::min(1.0, xi));
		double closest[3] = {
			piece.a[0] + xi * ab[0],
			piece.a[1] + xi * ab[1],
			piece.a[2] + xi * ab[2]};
		double delta[3] = {p[0] - closest[0], p[1] - closest[1], p[2] - closest[2]};
		const double distance = norm3(delta);
		if (distance < best_distance) {
			best_distance = distance;
			best = n;
			best_xi = xi;
		}
	}
	if (best < 0 || best_distance > 1e-4) {
		return false;
	}
	const Piece& piece = pieces_[best];
	const WireExport& view = exports_[piece.wire];
	int local = -1;
	for (int s = 0; s < static_cast<int>(view.piece_ids.size()); ++s) {
		if (view.piece_ids[s] == best) {
			local = s;
			break;
		}
	}
	if (local < 0 || !view.space) {
		return false;
	}
	mfem::IntegrationPoint ip;
	ip.Set1w(best_xi, 0.0);
	mfem::Vector shape(order_ + 1);
	view.space->GetFE(local)->CalcShape(ip, shape);
	const double* line_current = current_.HostRead();
	const double* line_charge = charge_.HostRead();
	current = 0.0;
	charge = 0.0;
	for (int i = 0; i < shape.Size(); ++i) {
		current += line_current[piece.offset + i] * shape(i);
		charge += line_charge[piece.offset + i] * shape(i);
	}
	return true;
}

void ThinWireCoupler::buildExportMeshes()
{
	const int nwire = static_cast<int>(model_.wireCurves().size());
	exports_.clear();
	exports_.resize(nwire);
	for (int p = 0; p < static_cast<int>(pieces_.size()); ++p) {
		exports_[pieces_[p].wire].piece_ids.push_back(p);
	}
	const int nloc = order_ + 1;
	for (int w = 0; w < nwire; ++w) {
		WireExport& view = exports_[w];
		const int nseg = static_cast<int>(view.piece_ids.size());
		if (nseg == 0) {
			continue;
		}
		view.mesh = std::make_unique<mfem::Mesh>(1, nseg + 1, nseg, 0, 3);
		for (int s = 0; s < nseg; ++s) {
			const Piece& piece = pieces_[view.piece_ids[s]];
			if (s == 0) {
				view.mesh->AddVertex(piece.a[0], piece.a[1], piece.a[2]);
			}
			view.mesh->AddVertex(piece.b[0], piece.b[1], piece.b[2]);
			view.mesh->AddSegment(s, s + 1, 1);
		}
		view.mesh->FinalizeTopology();
		view.collection = std::make_unique<mfem::DG_FECollection>(order_, 1, basis_type_);
		view.space = std::make_unique<mfem::FiniteElementSpace>(view.mesh.get(), view.collection.get());
		view.current = std::make_unique<mfem::GridFunction>(view.space.get());
		view.charge = std::make_unique<mfem::GridFunction>(view.space.get());
		*view.current = 0.0;
		*view.charge = 0.0;
		if (view.space->GetNDofs() != nseg * nloc) {
			throw std::runtime_error("Wire export space does not match the line degrees of freedom.");
		}
	}
}

void ThinWireCoupler::syncExportFields() const
{
	const int nloc = order_ + 1;
	const double* line_current = current_.HostRead();
	const double* line_charge = charge_.HostRead();
	for (WireExport& view : exports_) {
		if (!view.current) {
			continue;
		}
		for (int s = 0; s < static_cast<int>(view.piece_ids.size()); ++s) {
			const Piece& piece = pieces_[view.piece_ids[s]];
			for (int i = 0; i < nloc; ++i) {
				(*view.current)(s * nloc + i) = line_current[piece.offset + i];
				(*view.charge)(s * nloc + i) = line_charge[piece.offset + i];
			}
		}
	}
}

void ThinWireCoupler::writeParaView(const std::string& directory, int cycle, double time) const
{
	if (mfem::Mpi::WorldRank() != 0 || exports_.empty()) {
		return;
	}
	const std::filesystem::path volume_dir(directory);
	const std::string volume_name = volume_dir.filename().string();
	std::string stem = volume_name;
	const std::string suffix = ".msh";
	if (stem.size() > suffix.size() && stem.compare(stem.size() - suffix.size(), suffix.size(), suffix) == 0) {
		stem.erase(stem.size() - suffix.size());
	}
	const std::filesystem::path prefix = volume_dir;

	syncExportFields();
	for (int w = 0; w < static_cast<int>(exports_.size()); ++w) {
		WireExport& view = exports_[w];
		if (!view.mesh) {
			continue;
		}
		if (!view.paraview) {
			std::ostringstream index;
			index << std::setw(2) << std::setfill('0') << w;
			const std::string name = stem + "_wire" + index.str() + ".msh";
			view.paraview = std::make_unique<mfem::ParaViewDataCollection>(name, view.mesh.get());
			view.paraview->SetPrefixPath(prefix.string());
			view.paraview->RegisterField("current", view.current.get());
			view.paraview->RegisterField("charge", view.charge.get());
			view.paraview->SetHighOrderOutput(order_ > 1);
			view.paraview->SetLevelsOfDetail(std::max(order_, 1));
			view.paraview->SetDataFormat(mfem::VTKFormat::BINARY);
		}
		view.paraview->SetCycle(cycle);
		view.paraview->SetTime(time);
		view.paraview->Save();
	}
}

void ThinWireCoupler::packDevice()
{
	const int npieces = static_cast<int>(pieces_.size());
	const int nloc = order_ + 1;
	if (npieces == 0) {
		return;
	}
	if (nloc > kMaxLineDofs) {
		throw std::runtime_error(
			"Thin-wire CUDA path supports polynomial order at most 15.");
	}
	nblocks_ = 0;
	nvol_ = 0;
	for (const Piece& piece : pieces_) {
		for (const TubeBlock& block : piece.blocks) {
			if (nvol_ == 0) {
				nvol_ = block.vol_dofs.Size();
			}
			if (block.vol_dofs.Size() != nvol_ || nvol_ <= 0) {
				throw std::runtime_error(
					"Thin-wire CUDA path expects equal-order tetrahedra in every tube.");
			}
			++nblocks_;
		}
	}

	block_piece_.SetSize(std::max(nblocks_, 1));
	block_piece_.UseDevice(true);
	vol_dofs_.SetSize(std::max(nblocks_, 1) * std::max(nvol_, 1));
	vol_dofs_.UseDevice(true);
	int* block_piece = block_piece_.HostWrite();
	int* dofs = vol_dofs_.HostWrite();
	const int es_stride = nloc * std::max(nvol_, 1);
	es_.SetSize(std::max(nblocks_, 1) * 3 * es_stride);
	es_.UseDevice(true);
	double* es = es_.HostWrite();
	const int j_stride = std::max(nvol_, 1) * nloc;
	j_mat_.SetSize(std::max(nblocks_, 1) * 3 * j_stride);
	j_mat_.UseDevice(true);
	double* jmat = j_mat_.HostWrite();
	int block_index = 0;
	for (int p = 0; p < npieces; ++p) {
		for (const TubeBlock& block : pieces_[static_cast<std::size_t>(p)].blocks) {
			block_piece[block_index] = p;
			for (int i = 0; i < nvol_; ++i) {
				dofs[block_index * nvol_ + i] = block.vol_dofs[i];
			}
			const mfem::DenseMatrix* es_blocks[3] = {&block.es_x, &block.es_y, &block.es_z};
			const mfem::DenseMatrix* j_blocks[3] = {&block.j_x, &block.j_y, &block.j_z};
			for (int d = 0; d < 3; ++d) {
				for (int row = 0; row < nloc; ++row) {
					for (int col = 0; col < nvol_; ++col) {
						es[(block_index * 3 + d) * es_stride + row * nvol_ + col] =
							(*es_blocks[d])(row, col);
					}
				}
				for (int row = 0; row < nvol_; ++row) {
					for (int col = 0; col < nloc; ++col) {
						jmat[(block_index * 3 + d) * j_stride + row * nloc + col] =
							(*j_blocks[d])(row, col);
					}
				}
			}
			++block_index;
		}
	}

	mass_inv_.SetSize(npieces * nloc * nloc);
	mass_inv_.UseDevice(true);
	double* minv = mass_inv_.HostWrite();
	stiffness_.SetSize(npieces * nloc * nloc);
	stiffness_.UseDevice(true);
	double* stiff = stiffness_.HostWrite();
	for (int p = 0; p < npieces; ++p) {
		const Piece& piece = pieces_[static_cast<std::size_t>(p)];
		for (int row = 0; row < nloc; ++row) {
			for (int col = 0; col < nloc; ++col) {
				const int slot = (p * nloc + row) * nloc + col;
				minv[slot] = piece.mass_inv(row, col);
				stiff[slot] = piece.stiffness(row, col);
			}
		}
	}

	shape_left_.SetSize(nloc);
	shape_right_.SetSize(nloc);
	shape_left_.UseDevice(true);
	shape_right_.UseDevice(true);
	double* left = shape_left_.HostWrite();
	double* right = shape_right_.HostWrite();
	for (int i = 0; i < nloc; ++i) {
		left[i] = pieces_[0].shape_left(i);
		right[i] = pieces_[0].shape_right(i);
	}

	v2_.SetSize(npieces);
	inductance_.SetSize(npieces);
	resistance_.SetSize(npieces);
	v2_.UseDevice(true);
	inductance_.UseDevice(true);
	resistance_.UseDevice(true);
	double* v2 = v2_.HostWrite();
	double* inductance = inductance_.HostWrite();
	double* resistance = resistance_.HostWrite();

	offset_.SetSize(npieces);
	left_neighbor_.SetSize(npieces);
	right_neighbor_.SetSize(npieces);
	left_pec_.SetSize(npieces);
	right_pec_.SetSize(npieces);
	offset_.UseDevice(true);
	left_neighbor_.UseDevice(true);
	right_neighbor_.UseDevice(true);
	left_pec_.UseDevice(true);
	right_pec_.UseDevice(true);
	int* offset = offset_.HostWrite();
	int* left_neighbor = left_neighbor_.HostWrite();
	int* right_neighbor = right_neighbor_.HostWrite();
	int* left_pec = left_pec_.HostWrite();
	int* right_pec = right_pec_.HostWrite();
	for (int p = 0; p < npieces; ++p) {
		v2[p] = pieces_[p].v2;
		inductance[p] = pieces_[p].inductance;
		resistance[p] = pieces_[p].resistance;
		offset[p] = pieces_[p].offset;
		left_neighbor[p] = pieces_[p].left_neighbor;
		right_neighbor[p] = pieces_[p].right_neighbor;
		left_pec[p] = pieces_[p].left_pec ? 1 : 0;
		right_pec[p] = pieces_[p].right_pec ? 1 : 0;
	}

	drive_.SetSize(current_.Size());
	drive_.UseDevice(true);
}

}  // namespace maxwell
