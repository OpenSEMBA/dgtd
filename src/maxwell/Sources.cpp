#include "Sources.h"

#include "PhysicalConstants.h"
#include "Calculus.h"

namespace maxwell {

constexpr double TOLERANCE = 10.0*DBL_EPSILON;

mfem::DenseMatrix getRotationMatrix(const Source::CartesianAngles& angles_)
{
	std::array<mfem::DenseMatrix, 3> rotMat;
	double rotationAxisX[3] = { 1.0,0.0,0.0 }, rotationAxisY[3] = { 0.0,1.0,0.0 }, rotationAxisZ[3] = { 0.0,0.0,1.0 };
	mfem::NURBSPatch::Get3DRotationMatrix(rotationAxisX, angles_[X], 1.0, rotMat[X]);
	mfem::NURBSPatch::Get3DRotationMatrix(rotationAxisY, angles_[Y], 1.0, rotMat[Y]);
	mfem::NURBSPatch::Get3DRotationMatrix(rotationAxisZ, angles_[Z], 1.0, rotMat[Z]);
	mfem::DenseMatrix tMat(3), res(3);
	mfem::Mult(rotMat[X], rotMat[Y], tMat);
	mfem::Mult(tMat     , rotMat[Z], res);
	return res;
}

mfem::Vector rotateAroundAxis(
	const mfem::Vector& v, const Source::CartesianAngles& angles_)
{
	auto rotMat{ getRotationMatrix(angles_) };

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
	magnitude_{ f.clone() },
	fieldType_{ fT },
	polarization_{ p },
	center_{ centerIn },
	angles_{ anglesIn }
{
	assert(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE);
}

InitialField::InitialField(const InitialField& rhs) :
	magnitude_{ rhs.magnitude_->clone() },
	fieldType_{ rhs.fieldType_ },
	polarization_{ rhs.polarization_ },
	center_{ rhs.center_ },
	angles_{ rhs.angles_ }
{}

std::unique_ptr<Source> InitialField::clone() const
{
	return std::make_unique<InitialField>(*this);
}

double InitialField::eval(
	const Position& p, const Time& t,
	const FieldType& f, const Direction& d) const
{
	assert(p.Size() == center_.Size());
	
	if (f != fieldType_) {
		return 0.0;
	}

	Position pos(p.Size());
	for (int i{ 0 }; i < p.Size(); ++i) {
		pos[i] = p[i] - center_[i];
	}
	if (angles_[0] != 0.0 || angles_[1] != 0.0 || angles_[2] != 0.0) {
		pos = rotateAroundAxis(pos, angles_);
	}
	
	return magnitude_->eval(pos) * polarization_[d];
}

Planewave::Planewave(
	const MathFunction& mag, 
	const Polarization& p, 
	const Propagation& dir):
	magnitude_{ mag.clone() },
	polarization_{ p },
	propagation_{ dir }
{}

Planewave::Planewave(const Planewave& rhs) :
	magnitude_{ rhs.magnitude_->clone() },
	polarization_{ rhs.polarization_ },
	propagation_{ rhs.propagation_ }
{
	assert(polarization_.Norml2() <= TOLERANCE);
	assert(propagation_.Norml2() <= TOLERANCE);
}

std::unique_ptr<Source> Planewave::clone() const
{
	return std::make_unique<Planewave>(*this);
}

double Planewave::eval(
	const Position& p, const Time& t,
	const FieldType& f, const Direction& d) const
{
	assert(f == E || f == H);
	assert(d == X || d == Y || d == Z);
	
	mfem::Vector fieldPol(3);
	if (f == E) {
		fieldPol = polarization_;
	}
	else {
		fieldPol = crossProduct(propagation_, polarization_);
	}

	mfem::Vector delayedPosition(
		{ p[d] - t * propagation_[d] / physicalConstants::speedOfLight }
	);
	return magnitude_->eval(delayedPosition) * fieldPol[d];
}

}