#pragma once

#include <mfem.hpp>

namespace maxwell {

static mfem::Vector crossProduct(const mfem::Vector& va, const mfem::Vector& vb)
{
	assert(va.Size() == 3);
	assert(vb.Size() == 3);
	mfem::Vector V(3);
	V[0] = va[1] * vb[2] - va[2] * vb[1];
	V[1] = va[2] * vb[0] - va[0] * vb[2];
	V[2] = va[0] * vb[1] - va[1] * vb[0];
	return V;
}

static mfem::Vector unitVec(int d)
{
	
	mfem::Vector r(3);
	r = 0.0;
	r[d] = 1.0;
	return r;
}

static mfem::Vector minusUnitVec(int d)
{
	auto res{unitVec(d)};
	res *= -1.0;
	return res;
}




}