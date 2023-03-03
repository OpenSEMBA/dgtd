#include "Sources.h"

namespace maxwell {

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
	//if p.size() == 2 then apply 2D rotMat, else apply 3D (on Z?) Incomplete TODO
	assert(p.Size() == center.Size());
	Position pos(p.Size());
	for (int i{ 0 }; i < p.Size(); ++i) {
		pos[i] = p[i] - center[i];
	}
	if (rotAngle != 0) {
		auto tPos = pos;
		mfem::DenseMatrix rotMat;
		double rotationAxis[3] = { 0.0, 0.0, 1.0 };
		mfem::NURBSPatch::Get3DRotationMatrix(rotationAxis, -M_PI_4, 1.0, rotMat);
		mfem::Vector pos3D(3), newPos3D(3);
		pos3D[0] = pos[0];
		pos3D[1] = pos[1];
		pos3D[2] = 0.0;
		rotMat.Mult(pos3D, newPos3D);
		pos[0] = newPos3D[0];
		pos[1] = newPos3D[1];
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