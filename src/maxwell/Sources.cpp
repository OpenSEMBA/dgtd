#include "Sources.h"

namespace maxwell {

mfem::Vector rotateAroundZAxis(
	const mfem::Vector& v, const double& angleRads)
{
	mfem::DenseMatrix rotMat;
	double rotationAxis[3] = { 0.0, 0.0, 1.0 };
	mfem::NURBSPatch::Get3DRotationMatrix(rotationAxis, angleRads, 1.0, rotMat);

	mfem::Vector pos3D(3), newPos3D(3);
	pos3D = 0.0;
	for (auto d{ 0 }; d < 3; d++) {
		pos3D[d] = v[d];
	}
	
	rotMat.Mult(pos3D, newPos3D);

	mfem::Vector res(v.Size());
	for (auto d{ 0 }; d < v.Size(); d++) {
		res[d] = v[d];
	}
	return res;
}

InitialField::InitialField(
	const MathFunction& f, 
	const FieldType& fT, 
	const Polarization& p,
	const Position& centerIn,
	const double rotAngleIn ) :
	function_{ f.clone() },
	fieldType{ fT },
	polarization{ p },
	center{ centerIn },
	rotAngle{ rotAngleIn }
{}

InitialField::InitialField(const InitialField& rhs) :
	function_{ rhs.function_->clone() },
	fieldType{ rhs.fieldType },
	polarization{ rhs.polarization },
	center{ rhs.center },
	rotAngle{ rhs.rotAngle }
{}

std::unique_ptr<Source> InitialField::clone() const
{
	return std::make_unique<InitialField>(*this);
}

double InitialField::eval(const Position& p, Time t) const
{
	assert(p.Size() == center.Size());
	Position pos(p.Size());
	for (int i{ 0 }; i < p.Size(); ++i) {
		pos[i] = p[i] - center[i];
	}
	if (rotAngle != 0) {
		pos = rotateAroundZAxis(pos, -M_PI_4);
	}
	return function_->eval(pos, t);
}

PlaneWave::PlaneWave(const MathFunction& f, const Polarization& p):
	function_{ f.clone() },
	polarization{ p }
{}

PlaneWave::PlaneWave(const PlaneWave& rhs) :
	function_{ rhs.function_->clone() },
	polarization{ rhs.polarization }
{}

std::unique_ptr<Source> PlaneWave::clone() const
{
	return std::make_unique<PlaneWave>(*this);
}

double PlaneWave::eval(const Position& p, Time t) const
{
	return function_->eval(p, t);
}

}