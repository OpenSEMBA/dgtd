#include "Sources.h"

namespace maxwell {

mfem::DenseMatrix getRotationMatrix(const Source::CartesianAngles& angles)
{
	std::array<mfem::DenseMatrix, 3> rotMat;
	double rotationAxisX[3] = { 1.0,0.0,0.0 }, rotationAxisY[3] = { 0.0,1.0,0.0 }, rotationAxisZ[3] = { 0.0,0.0,1.0 };
	mfem::NURBSPatch::Get3DRotationMatrix(rotationAxisX, angles[X], 1.0, rotMat[X]);
	mfem::NURBSPatch::Get3DRotationMatrix(rotationAxisY, angles[Y], 1.0, rotMat[Y]);
	mfem::NURBSPatch::Get3DRotationMatrix(rotationAxisZ, angles[Z], 1.0, rotMat[Z]);
	mfem::DenseMatrix tMat(3), res(3);
	mfem::Mult(rotMat[X], rotMat[Y], tMat);
	mfem::Mult(tMat     , rotMat[Z], res);
	return res;
}

mfem::Vector rotateAroundAxis(
	const mfem::Vector& v, const Source::CartesianAngles& angles)
{
	auto rotMat{ getRotationMatrix(angles) };

	mfem::Vector pos3D(3), newPos3D(3);
	pos3D = 0.0;
	for (auto d{ 0 }; d < v.Size(); d++) {
		pos3D[d] = v[d];
	}

	rotMat.Mult(pos3D, newPos3D);

	mfem::Vector res(v.Size());
	for (auto d{ 0 }; d < v.Size(); d++) {
		res[d] = newPos3D[d];
	}
	return res;
}

InitialField::InitialField(
	const MathFunction& f, 
	const FieldType& fT, 
	const Polarization& p,
	const Position& centerIn,
	const CartesianAngles anglesIn ) :
	function_{ f.clone() },
	fieldType{ fT },
	polarization{ p },
	center{ centerIn },
	angles{ anglesIn }
{}

InitialField::InitialField(const InitialField& rhs) :
	function_{ rhs.function_->clone() },
	fieldType{ rhs.fieldType },
	polarization{ rhs.polarization },
	center{ rhs.center },
	angles{ rhs.angles }
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
	if (angles[0] != 0.0 || angles[1] != 0.0 || angles[2] != 0.0) {
		auto tPos = rotateAroundAxis(pos, angles);
		return function_->eval(tPos, t);
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