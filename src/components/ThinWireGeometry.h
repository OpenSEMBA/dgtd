#pragma once

#include <array>
#include <vector>

namespace maxwell {

struct WireSpec {
	std::vector<int> tags;
	double radius = 0.0;
	double resistance_per_meter = 0.0;
	double inductance_per_meter = 0.0;
	// Transverse length Δ. The tube radius is 1.7Δ. Zero uses the host tet mean edge.
	double transverse_step = 0.0;
};

struct WirePolyline {
	std::vector<int> tags;
	std::vector<std::array<double, 3>> points;
	bool closed = false;
};

// One tetrahedron-sized piece of a thin wire. The chain is identical on every rank.
// `host` is the serial-mesh element index.
struct WirePieceGeom {
	int host = -1;
	int wire = 0;
	int left_neighbor = -1;
	int right_neighbor = -1;
	bool left_pec = false;
	bool right_pec = false;
	double a[3] = {0.0, 0.0, 0.0};
	double b[3] = {0.0, 0.0, 0.0};
	double length = 0.0;
	double tangent[3] = {0.0, 0.0, 0.0};
	double inductance = 1.0;
	double resistance = 0.0;
	double v2 = 1.0;
	double radius = 0.0;
	double tube_radius = 0.0;
};

}
