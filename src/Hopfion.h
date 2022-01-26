#pragma once

#include <array>
#include <complex>
#include <stdlib.h>
#include <vector>
#include <mfem.hpp>


class Hopfion {
public:

	typedef std::array<double, 3> Vec3;
	typedef std::pair<Vec3, Vec3> FieldEH;
	
	Hopfion(std::size_t p, std::size_t q);

	FieldEH evaluate(double time, Vec3 position) const;


private:

	std::size_t p, q;


    double FieldEx(double x, double y, double z, double t) const{
        double Ex;
        std::complex<double> I(0, 1);
        std::complex<double> Fx;

        Fx = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) * (-4 * x * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2) + 2.0 * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * x / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2)) - 2 * x * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2));
        Ex = real(Fx);

        return Ex;
    }

    double FieldHx(double x, double y, double z, double t) const {
        double Hx;
        std::complex<double> I(0, 1);
        std::complex<double> Fx;

        Fx = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) * (-4 * x * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2) + 2.0 * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * x / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2)) - 2 * x * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2));
        Hx = imag(Fx);

        return Hx;
    }


    double FieldEy(double x, double y, double z, double t) const {
        double Ey;
        std::complex<double> I(0, 1);
        std::complex<double> Fy;

        Fy = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) * (-4 * y * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 2.0 * I * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * y / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2 * y * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Ey = real(Fy);

        return Ey;
    }

    double FieldHy(double x, double y, double z, double t) const {
        double Hy;
        std::complex<double> I(0, 1);
        std::complex<double> Fy;

        Fy = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) * (-4 * y * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 2.0 * I * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * y / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2 * y * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Hy = imag(Fy);

        return Hy;
    }


    double FieldEz(double x, double y, double z, double t) const {
        double Ez;
        std::complex<double> I(0, 1);
        std::complex<double> Fz;

        Fz = I * (-4 * z * (x - I * y) * (-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 4.0 * (x - I * y) * (-I + t) * ((2.0 * I + 2.0 * z) / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2.0 * z * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Ez = real(Fz);

        return Ez;
    }

    double FieldHz(double x, double y, double z, double t) const {
        double Hz;
        std::complex<double> I(0, 1);
        std::complex<double> Fz;

        Fz = I * (-4 * z * (x - I * y) * (-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 4.0 * (x - I * y) * (-I + t) * ((2.0 * I + 2.0 * z) / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2.0 * z * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Hz = imag(Fz);

        return Hz;
    }

};


std::vector<int> mapQuadElementTopLeftVertex(const mfem::Mesh& mesh);