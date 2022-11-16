#pragma once

#include <mfem.hpp>

namespace AnalyticalFunctions3D {
	mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;

	double gaussianFunction(const mfem::Vector& pos)
	{
		mfem::Vector normalizedPos(3);
		for (size_t i = 0; i < 3; i++) {
			double center = (meshBoundingBoxMin[i] + meshBoundingBoxMax[i]) * 0.5;
			normalizedPos[i] = 2 * (pos[i] - center) / (meshBoundingBoxMax[i] - meshBoundingBoxMin[i]);
		}

		return exp(-20. * (pow(normalizedPos[0], 2) + pow(normalizedPos[1], 2) + pow(normalizedPos[2], 2)));
	}
}