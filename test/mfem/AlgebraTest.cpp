#include <gtest/gtest.h>

#include <iostream>
#include <filesystem>
#include <mfem.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include "math/PhysicalConstants.h"

#include <math.h>
#include <complex>
#include <fftw3.h>

using namespace mfem;
class AlgebraTest : public ::testing::Test 
{
public:

	double func_exp_real_part_2D(const Vector& x, const double freq, const double phi)
	{
		//angulo viene dado por x[0], x[1] y 0.0, 0.0. No es el angulo donde observo, es el angulo que forma el punto y el angulo de observacion en un sistema centrado en el punto.
		auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) 
			/ sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
		return cos(2.0 * M_PI * (maxwell::physicalConstants::speedOfLight_SI / freq) * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
	}

	double func_exp_imag_part_2D(const Vector& x, const double freq, const double phi)
	{
		auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
		return sin(2.0 * M_PI * (maxwell::physicalConstants::speedOfLight_SI / freq) * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
	}

	std::complex<double> func_exp_2D(const Vector& x, const double freq, const double phi)
	{
		auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
		return exp(std::complex(0.0, 2.0 * M_PI * (maxwell::physicalConstants::speedOfLight_SI / freq) * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle)));
	}
};

std::complex<double> operator"" _i(long double x)
{
	return std::complex<double>(0.0, x);
}

TEST_F(AlgebraTest, calcRealEigenvalues)
{

	Eigen::Matrix3d matrix{
		{1.0, 0.0, 0.0},
		{0.0, 2.0, 0.0},
		{0.0, 0.0,-3.0}
	};

	Eigen::Vector<double, 3> expectedEVs{
		{1.0, 2.0, -3.0}
	};

	EXPECT_EQ(matrix.eigenvalues(), expectedEVs);
}

TEST_F(AlgebraTest, calcComplexEigenvalues)
{
	Eigen::Matrix2d matrix{
		{  3, -2},
		{  4, -1}
	};

	Eigen::Vector<std::complex<double>, 2> expectedEVs{
		{1.0 + 2.0_i, 1.0 - 2.0_i}
	};

	EXPECT_EQ(matrix.eigenvalues(), expectedEVs);

}

TEST_F(AlgebraTest, checkSparseMatrixVectorProduct)
{
	int a{ 4 }, b{ 4 };
	Eigen::SparseMatrix<double> sparse(a,b);
	sparse.insert(0, 0) = 1.0;
	sparse.insert(1, 0) = 2.0;
	sparse.insert(0, 2) = 3.0;
	sparse.insert(2, 3) = 4.0;
	sparse.makeCompressed();
	int c = 0;
	Eigen::Vector4d vec({ 1.0,1.0,1.0,1.0 });
	
	EXPECT_NO_THROW( sparse * vec );
}

TEST_F(AlgebraTest, checkEigenMatrixVectorProduct)
{
	int a{ 2 }, b{ 2 };
	Eigen::SparseMatrix<double> sparse(a, b);
	sparse.insert(0, 0) = 2;
	sparse.insert(0, 1) = -12;
	sparse.insert(1, 0) = 1;
	sparse.insert(1, 1) = -5;
	Eigen::VectorXd it(2);
	it[0] = 1;
	it[1] = 1;

	EXPECT_NO_THROW( sparse * it);
}

TEST_F(AlgebraTest, custom_rcs_exponential_function)
{
	Vector vec({ 1.5, 1.5 });
	double tol{ 1e-5 };

	EXPECT_NEAR(func_exp_real_part_2D(vec, 100000.0, 0.0),  0.757829, tol);
	EXPECT_NEAR(func_exp_imag_part_2D(vec, 100000.0, 0.0), -0.652453, tol);
	EXPECT_NEAR(func_exp_real_part_2D(vec, 100000.0, M_PI), 0.757829, tol);
	EXPECT_NEAR(func_exp_imag_part_2D(vec, 100000.0, M_PI), 0.652453, tol);

	EXPECT_NEAR(func_exp_2D(vec, 100000.0, 0.0).real(), func_exp_real_part_2D(vec, 100000.0, 0.0), tol);
	EXPECT_NEAR(func_exp_2D(vec, 100000.0, 0.0).imag(), func_exp_imag_part_2D(vec, 100000.0, 0.0), tol);
	EXPECT_NEAR(func_exp_2D(vec, 100000.0, M_PI).real(), func_exp_real_part_2D(vec, 100000.0, M_PI), tol);
	EXPECT_NEAR(func_exp_2D(vec, 100000.0, M_PI).imag(), func_exp_imag_part_2D(vec, 100000.0, M_PI), tol);
}

TEST_F(AlgebraTest, custom_dft_method)
{
	std::vector<std::vector<double>> time_fields{ { 2.5, 3.5 }, { 2.0, 4.0 }, { 1.5, 4.5 } };
	std::vector<double> times{ {0.0, 5e-3, 1e-2} };
	std::vector<double> frequencies{ 300e6, 500e6, 700e6, 900e6 };
	std::vector<std::vector<std::complex<double>>> freq_fields(frequencies.size());

	for (int t{ 0 }; t < times.size(); ++t) {
		for (int f{ 0 }; f < frequencies.size(); ++f) {
			freq_fields[f].resize(time_fields[t].size());
			auto time_const{ std::exp(std::complex<double>(0.0, -2.0 * M_PI * frequencies[f] * (times[t] / maxwell::physicalConstants::speedOfLight_SI))) };
			for (int v{ 0 }; v < time_fields[t].size(); ++v) {
				freq_fields[f][v] += std::complex<double>(time_fields[t][v], 0.0) * time_const;
			}
		}
	}

	double tol{ 1e-4 };
	std::complex<double> expected_300_0(5.99605, -0.157116);
	std::complex<double> expected_300_1(11.9891, -0.408483);
	std::complex<double> expected_500_0(5.98903, -0.261645);
	std::complex<double> expected_500_1(11.9698, -0.680191);
	std::complex<double> expected_700_0(5.97851, -0.365853);
	std::complex<double> expected_700_1(11.9409, -0.950981);
	std::complex<double> expected_900_0(5.96451, -0.469611);
	std::complex<double> expected_900_1(11.9024, -1.22049);

	EXPECT_NEAR(freq_fields[0][0].real(), expected_300_0.real(), tol);
	EXPECT_NEAR(freq_fields[0][0].imag(), expected_300_0.imag(), tol);
	EXPECT_NEAR(freq_fields[0][1].real(), expected_300_1.real(), tol);
	EXPECT_NEAR(freq_fields[0][1].imag(), expected_300_1.imag(), tol);
	EXPECT_NEAR(freq_fields[1][0].real(), expected_500_0.real(), tol);
	EXPECT_NEAR(freq_fields[1][0].imag(), expected_500_0.imag(), tol);
	EXPECT_NEAR(freq_fields[1][1].real(), expected_500_1.real(), tol);
	EXPECT_NEAR(freq_fields[1][1].imag(), expected_500_1.imag(), tol);
	EXPECT_NEAR(freq_fields[2][0].real(), expected_700_0.real(), tol);
	EXPECT_NEAR(freq_fields[2][0].imag(), expected_700_0.imag(), tol);
	EXPECT_NEAR(freq_fields[2][1].real(), expected_700_1.real(), tol);
	EXPECT_NEAR(freq_fields[2][1].imag(), expected_700_1.imag(), tol);
	EXPECT_NEAR(freq_fields[3][0].real(), expected_900_0.real(), tol);
	EXPECT_NEAR(freq_fields[3][0].imag(), expected_900_0.imag(), tol);
	EXPECT_NEAR(freq_fields[3][1].real(), expected_900_1.real(), tol);
	EXPECT_NEAR(freq_fields[3][1].imag(), expected_900_1.imag(), tol);

}

const double getTime(const std::string& timePath)
{
	std::ifstream timeFile(timePath);
	if (!timeFile) {
		throw std::runtime_error("File could not be opened in getTime.");
	}
	std::string timeString;
	std::getline(timeFile, timeString);
	return std::stod(timeString);
}

static void performDFT(const std::map<double, std::vector<double>>& dataMap,
	const std::vector<double>& times,
	std::vector<double>& frequencies,
	std::vector<std::vector<double>>& amplitudeSpectra) {
	int numTimePoints = int(times.size());
	int numFrequencyPoints = numTimePoints / 2 + 1;

	// Determine the time step (assuming equally spaced samples)
	double timeStep = times[1] - times[0];

	// Determine the maximum frequency based on the sampling rate (Nyquist frequency)
	double maxFrequency = 1.0 / (2 * timeStep);

	// Interpolate data onto a regular grid in the time domain
	std::vector<double> interpolatedData(numTimePoints);
	for (int i = 0; i < numTimePoints; ++i) {
		interpolatedData[i] = dataMap.lower_bound(times[i])->second[0]; // Assuming one component for simplicity
	}

	// Prepare input and output arrays for FFTW
	fftw_complex* out = (fftw_complex*)fftw_malloc(numFrequencyPoints * sizeof(fftw_complex));
	double* in = (double*)fftw_malloc(numTimePoints * sizeof(double));

	// Create FFTW plan for forward transform
	fftw_plan plan = fftw_plan_dft_r2c_1d(numTimePoints, in, out, FFTW_ESTIMATE);

	// Populate input array with interpolated data
	for (int i = 0; i < numTimePoints; ++i) {
		in[i] = interpolatedData[i];
	}

	// Perform forward DFT
	fftw_execute(plan);

	// Extract amplitudes from the frequency domain representation
	amplitudeSpectra.clear();
	for (int i = 0; i < numFrequencyPoints; ++i) {
		std::vector<double> amplitudes(dataMap.begin()->second.size(), 0.0); // Initialize with zeros
		amplitudes[0] = out[i][0]; // Real part of the spectrum
		if (i > 0 && i < numFrequencyPoints - 1) {
			amplitudes[i] = out[i][1]; // Imaginary part of the spectrum (except at Nyquist frequency)
			amplitudes[numTimePoints - i] = out[i][1]; // Imaginary part at mirrored frequency
		}
		else if (i == numFrequencyPoints - 1) {
			amplitudes[i] = out[i][1]; // Imaginary part at Nyquist frequency
		}
		amplitudeSpectra.push_back(amplitudes);
	}

	// Calculate frequencies
	frequencies.clear();
	for (int i = 0; i < numFrequencyPoints; ++i) {
		frequencies.push_back(i * maxFrequency / numTimePoints);
	}

	// Clean up
	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);
}