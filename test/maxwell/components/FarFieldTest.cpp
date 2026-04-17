#include "gtest/gtest.h"

#include "components/FarField.h"

#include <cmath>
#include <complex>
#include <map>
#include <vector>

using namespace maxwell;

class FarFieldTest : public ::testing::Test {};

// ----------------------------------------------------------------
// calcPsiAngle3D
// ----------------------------------------------------------------

TEST_F(FarFieldTest, calcPsiAngle3D_alignedWithNorthPole)
{
	// obs direction: theta=0, phi=0 -> (0, 0, r). p along +z -> psi = 0
	SphericalAngles angles{ 0.0, 0.0 };
	mfem::Vector p({ 0.0, 0.0, 1.0 });
	EXPECT_NEAR(0.0, calcPsiAngle3D(p, angles), 1e-10);
}

TEST_F(FarFieldTest, calcPsiAngle3D_perpendicularToNorthPole)
{
	// obs direction: theta=0, phi=0 -> (0, 0, r). p along +x -> psi = pi/2
	SphericalAngles angles{ 0.0, 0.0 };
	mfem::Vector p({ 1.0, 0.0, 0.0 });
	EXPECT_NEAR(M_PI / 2.0, calcPsiAngle3D(p, angles), 1e-10);
}

TEST_F(FarFieldTest, calcPsiAngle3D_equatorialObservation)
{
	// obs direction: theta=pi/2, phi=0 -> (r, 0, 0). p along +x -> psi = 0
	SphericalAngles angles{ M_PI / 2.0, 0.0 };
	mfem::Vector p({ 1.0, 0.0, 0.0 });
	EXPECT_NEAR(0.0, calcPsiAngle3D(p, angles), 1e-10);
}

// ----------------------------------------------------------------
// complexInnerProduct
// ----------------------------------------------------------------

TEST_F(FarFieldTest, complexInnerProduct_realVectors)
{
	// (2+0j)*(3+0j) = 6+0j
	ComplexVector a = { { 2.0, 0.0 } };
	ComplexVector b = { { 3.0, 0.0 } };
	auto res = complexInnerProduct(a, b);
	EXPECT_NEAR(6.0, res.real(), 1e-12);
	EXPECT_NEAR(0.0, res.imag(), 1e-12);
}

TEST_F(FarFieldTest, complexInnerProduct_complexMultiplication)
{
	// (1+2j)*(3+4j) = 3+4j+6j+8j^2 = -5+10j
	ComplexVector a = { { 1.0, 2.0 } };
	ComplexVector b = { { 3.0, 4.0 } };
	auto res = complexInnerProduct(a, b);
	EXPECT_NEAR(-5.0, res.real(), 1e-12);
	EXPECT_NEAR(10.0, res.imag(), 1e-12);
}

TEST_F(FarFieldTest, complexInnerProduct_sizeMismatchThrows)
{
	ComplexVector a = { { 1.0, 0.0 }, { 2.0, 0.0 } };
	ComplexVector b = { { 1.0, 0.0 } };
	EXPECT_THROW(complexInnerProduct(a, b), std::runtime_error);
}

// ----------------------------------------------------------------
// func_exp_real_part_3D / func_exp_imag_part_3D
// ----------------------------------------------------------------

TEST_F(FarFieldTest, funcExp3D_antiAligned)
{
	// obs along +z (theta=0), p along -z -> psi=pi -> cos(psi)=-1
	// freq=1/(2*pi) -> wavenumber=1, |p|=1 -> rad_term=-1
	// real=cos(-1)=cos(1), imag=sin(-1)=-sin(1)
	SphericalAngles angles{ 0.0, 0.0 };
	mfem::Vector p({ 0.0, 0.0, -1.0 });
	const double freq = 1.0 / (2.0 * M_PI);
	EXPECT_NEAR(std::cos(1.0),  func_exp_real_part_3D(p, freq, angles), 1e-10);
	EXPECT_NEAR(-std::sin(1.0), func_exp_imag_part_3D(p, freq, angles), 1e-10);
}

TEST_F(FarFieldTest, funcExp3D_quarterWavelength)
{
	// obs along +z, p along +z -> psi=0, freq=0.25
	// landa=4, wavenumber=pi/2, rad_term=pi/2
	// real=cos(pi/2)=0, imag=sin(pi/2)=1
	SphericalAngles angles{ 0.0, 0.0 };
	mfem::Vector p({ 0.0, 0.0, 1.0 });
	EXPECT_NEAR(0.0, func_exp_real_part_3D(p, 0.25, angles), 1e-10);
	EXPECT_NEAR(1.0, func_exp_imag_part_3D(p, 0.25, angles), 1e-10);
}

TEST_F(FarFieldTest, funcExp3D_perpendicularObsAndPoint)
{
	// obs along +z (theta=0), p along +x -> psi=pi/2 -> cos(psi)=0 -> rad_term=0
	// => real=cos(0)=1, imag=sin(0)=0 regardless of frequency
	SphericalAngles angles{ 0.0, 0.0 };
	mfem::Vector p({ 1.0, 0.0, 0.0 });
	EXPECT_NEAR(1.0, func_exp_real_part_3D(p, 0.5, angles), 1e-12);
	EXPECT_NEAR(0.0, func_exp_imag_part_3D(p, 0.5, angles), 1e-12);
}

TEST_F(FarFieldTest, funcExp3D_alignedObsAndPoint)
{
	// obs along +z (theta=0), p along +z -> psi=0 -> rad_term = wavenumber * |p|
	// freq=1/(2*pi), speedOfLight=1 -> landa=2*pi, wavenumber=1 -> rad_term=1
	SphericalAngles angles{ 0.0, 0.0 };
	mfem::Vector p({ 0.0, 0.0, 1.0 });
	const double freq = 1.0 / (2.0 * M_PI);
	EXPECT_NEAR(std::cos(1.0), func_exp_real_part_3D(p, freq, angles), 1e-10);
	EXPECT_NEAR(std::sin(1.0), func_exp_imag_part_3D(p, freq, angles), 1e-10);
}

// ----------------------------------------------------------------
// func_exp_real_part_2D / func_exp_imag_part_2D
// ----------------------------------------------------------------

TEST_F(FarFieldTest, funcExp2D_alignedObsAndPoint)
{
	// obs along +x (theta=pi/2, phi=0), p along +x in XY -> psi=0
	// freq=1/(2*pi) -> wavenumber=1, |p|=1 -> rad_term=1
	SphericalAngles angles{ M_PI / 2.0, 0.0 };
	mfem::Vector p({ 1.0, 0.0, 0.0 });
	const double freq = 1.0 / (2.0 * M_PI);
	EXPECT_NEAR(std::cos(1.0), func_exp_real_part_2D(p, freq, angles), 1e-10);
	EXPECT_NEAR(std::sin(1.0), func_exp_imag_part_2D(p, freq, angles), 1e-10);
}

TEST_F(FarFieldTest, funcExp2D_quarterPeriodPoint)
{
	// obs along +x, p along +x, freq=0.25 -> wavenumber=2*pi*0.25=pi/2
	// |p|=2 -> rad_term=pi -> cos(pi)=-1, sin(pi)=0
	SphericalAngles angles{ M_PI / 2.0, 0.0 };
	mfem::Vector p({ 2.0, 0.0, 0.0 });
	EXPECT_NEAR(-1.0, func_exp_real_part_2D(p, 0.25, angles), 1e-10);
	EXPECT_NEAR(0.0,  func_exp_imag_part_2D(p, 0.25, angles), 1e-10);
}

// ----------------------------------------------------------------
// calculateDFT
// ----------------------------------------------------------------

TEST_F(FarFieldTest, calculateDFT_zeroFreqAtZeroTime)
{
	// freq=0, t=0 -> arg=0, w=(1,0) -> result = gf value
	mfem::Vector gf({ 3.0 });
	std::vector<Frequency> freqs = { 0.0 };
	auto res = calculateDFT(gf, freqs, 0.0);
	ASSERT_EQ(1u, res.size());
	ASSERT_EQ(1u, res[0].size());
	EXPECT_NEAR(3.0, res[0][0].real(), 1e-12);
	EXPECT_NEAR(0.0, res[0][0].imag(), 1e-12);
}

TEST_F(FarFieldTest, calculateDFT_halfPeriod)
{
	// freq=0.25, t=2.0 -> arg=2*pi*0.25*2=pi -> w=(cos(pi),-sin(pi))=(-1,0)
	// gf={1.0} -> res[0][0] = -1 + 0j
	mfem::Vector gf({ 1.0 });
	std::vector<Frequency> freqs = { 0.25 };
	auto res = calculateDFT(gf, freqs, 2.0);
	ASSERT_EQ(1u, res.size());
	ASSERT_EQ(1u, res[0].size());
	EXPECT_NEAR(-1.0, res[0][0].real(), 1e-10);
	EXPECT_NEAR(0.0,  res[0][0].imag(), 1e-10);
}

TEST_F(FarFieldTest, calculateDFT_multiDOF)
{
	// freq=0, t=0 -> w=(1,0) -> result = gf values unchanged
	mfem::Vector gf({ 2.0, 5.0 });
	std::vector<Frequency> freqs = { 0.0 };
	auto res = calculateDFT(gf, freqs, 0.0);
	ASSERT_EQ(1u, res.size());
	ASSERT_EQ(2u, res[0].size());
	EXPECT_NEAR(2.0, res[0][0].real(), 1e-12);
	EXPECT_NEAR(5.0, res[0][1].real(), 1e-12);
}

// ----------------------------------------------------------------
// evaluateGaussianVector
// ----------------------------------------------------------------

TEST_F(FarFieldTest, evaluateGaussianVector_peakAtMean)
{
	// t = |mean| -> exp(0) = 1.0
	std::vector<Time> times = { 2.0 };
	auto res = evaluateGaussianVector(times, 1.0, 2.0);
	ASSERT_EQ(1u, res.size());
	EXPECT_NEAR(1.0, res[0], 1e-12);
}

TEST_F(FarFieldTest, evaluateGaussianVector_oneSigmaFromMean)
{
	// t = |mean| + spread -> exp(-1/2) ≈ 0.6065
	std::vector<Time> times = { 3.0 };
	auto res = evaluateGaussianVector(times, 1.0, 2.0);
	ASSERT_EQ(1u, res.size());
	EXPECT_NEAR(std::exp(-0.5), res[0], 1e-12);
}

// ----------------------------------------------------------------
// trimLowMagFreqs
// ----------------------------------------------------------------

TEST_F(FarFieldTest, trimLowMagFreqs_noTrimNeeded)
{
	std::map<double, std::complex<double>> map{ {1.0, {0.5, 0.0}}, {2.0, {0.5, 0.0}} };
	std::vector<Frequency> freqs = { 1.0, 2.0 };
	trimLowMagFreqs(map, freqs);
	EXPECT_EQ(2u, freqs.size());
}

TEST_F(FarFieldTest, trimLowMagFreqs_trimsFromFirstLowValue)
{
	// f=1.0 is above tol, f=2.0 is below -> trim f=2.0 and everything after
	std::map<double, std::complex<double>> map{ {1.0, {0.5, 0.0}}, {2.0, {0.001, 0.0}} };
	std::vector<Frequency> freqs = { 1.0, 2.0 };
	trimLowMagFreqs(map, freqs);
	ASSERT_EQ(1u, freqs.size());
	EXPECT_NEAR(1.0, freqs[0], 1e-12);
}

TEST_F(FarFieldTest, trimLowMagFreqs_trimsAll)
{
	std::map<double, std::complex<double>> map{ {1.0, {0.001, 0.0}}, {2.0, {0.001, 0.0}} };
	std::vector<Frequency> freqs = { 1.0, 2.0 };
	trimLowMagFreqs(map, freqs);
	EXPECT_TRUE(freqs.empty());
}

// ----------------------------------------------------------------
// getNearToFarFieldMarker
// ----------------------------------------------------------------

TEST_F(FarFieldTest, getNearToFarFieldMarker_correctSize)
{
	auto marker = getNearToFarFieldMarker(201);
	EXPECT_EQ(201, marker.Size());
}

TEST_F(FarFieldTest, getNearToFarFieldMarker_correctIndex)
{
	// BdrCond::NearToFarField = 201 -> index 200 should be 1, rest 0
	auto marker = getNearToFarFieldMarker(201);
	const int idx = static_cast<int>(BdrCond::NearToFarField) - 1;
	EXPECT_EQ(1, marker[idx]);
	EXPECT_EQ(0, marker[0]);
	EXPECT_EQ(0, marker[1]);
}

// ----------------------------------------------------------------
// initAngles2FreqValues
// ----------------------------------------------------------------

TEST_F(FarFieldTest, initAngles2FreqValues_singleAngleSingleFreq)
{
	std::vector<Frequency> freqs = { 1.0 };
	std::vector<SphericalAngles> angles = { {0.0, 0.0} };
	auto res = initAngles2FreqValues(freqs, angles);
	ASSERT_EQ(1u, res.size());
	EXPECT_NEAR(0.0, res.at({0.0, 0.0}).at(1.0), 1e-12);
}

TEST_F(FarFieldTest, initAngles2FreqValues_multipleAnglesAndFreqs)
{
	std::vector<Frequency> freqs = { 1.0, 2.0 };
	std::vector<SphericalAngles> angles = { {0.0, 0.0}, {M_PI / 2.0, 0.0} };
	auto res = initAngles2FreqValues(freqs, angles);
	EXPECT_EQ(2u, res.size());
	EXPECT_EQ(2u, res.at({0.0, 0.0}).size());
	EXPECT_EQ(2u, res.at({M_PI / 2.0, 0.0}).size());
}

// ----------------------------------------------------------------
// PlaneWaveData
// ----------------------------------------------------------------

TEST_F(FarFieldTest, planeWaveData_notModulatedWhenFreqZero)
{
	PlaneWaveData pwd(1.0, 2.0);
	EXPECT_FALSE(pwd.isModulated());
	EXPECT_NEAR(1.0, pwd.spread, 1e-12);
	EXPECT_NEAR(2.0, pwd.mean, 1e-12);
}

TEST_F(FarFieldTest, planeWaveData_modulatedWhenFreqNonzero)
{
	PlaneWaveData pwd(1.0, 2.0, 3.0);
	EXPECT_TRUE(pwd.isModulated());
	EXPECT_NEAR(3.0, pwd.frequency, 1e-12);
}

// ----------------------------------------------------------------
// FreqFields
// ----------------------------------------------------------------

TEST_F(FarFieldTest, freqFields_appendAndAccessEx)
{
	FreqFields ff(2);
	ComplexVector cv = { {1.0, 0.0}, {2.0, 0.0} };
	ff.append(cv, "/Ex.gf", 0);
	ASSERT_EQ(2u, ff.Ex[0].size());
	EXPECT_NEAR(1.0, ff.Ex[0][0].real(), 1e-12);
	EXPECT_NEAR(2.0, ff.Ex[0][1].real(), 1e-12);
}

TEST_F(FarFieldTest, freqFields_normaliseFields)
{
	FreqFields ff(1);
	ComplexVector cv = { {4.0, 0.0} };
	ff.append(cv, "/Ex.gf", 0);
	ff.append(cv, "/Ey.gf", 0);
	ff.append(cv, "/Ez.gf", 0);
	ff.append(cv, "/Hx.gf", 0);
	ff.append(cv, "/Hy.gf", 0);
	ff.append(cv, "/Hz.gf", 0);
	ff.normaliseFields(2.0);
	EXPECT_NEAR(2.0, ff.Ex[0][0].real(), 1e-12);
	EXPECT_NEAR(2.0, ff.Hz[0][0].real(), 1e-12);
}
