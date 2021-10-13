#pragma once

#include <array>
#include <complex>

class SquareResonantCavity2D {
public:

	typedef std::array<double, 3> Vec3;
	typedef std::pair<Vec3, Vec3> FieldEH;

	SquareResonantCavity2D();


	FieldEH TheoField(double time, Vec3 position) const;
	FieldEH CompuField(double time) const;


private:

	// Funciones a utilizar para calcular los campos: completar correctamente.


    double theofieldEx(double x, double y, double z, double t) const {

        double Ex = 0; 

        return Ex;
    }

	double theofieldHx(double x, double y, double z, double t) const {

		double Hx = 0;

		return Hx;

};


	double theofieldEy(double x, double y, double z, double t) const {

		double Ey = 0;

		return Ey;
	}

	double theofieldHy(double x, double y, double z, double t) const {

		double Hy = 0;

		return Hy;
	}


	double theofieldHz(double x, double y, double z, double t) const {

		double Hz = 0;

		return Hz;
	}

	double theofieldEz(double x, double y, double z, double t) const {

		double Ez = 0;

		return Ez;
	}