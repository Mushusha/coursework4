#include "Quad.h"

std::vector<std::complex<double>> infQuad::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
	 FF.resize(4);
	 FF[0] = (1 + ksi / (1 - ksi)) * (1 - eta) / 2; // Q0
	 FF[1] = (-ksi / (1 - ksi)) * (1 - eta) / 2; // C0
	 FF[2] = (-ksi / (1 - ksi)) * (1 + eta) / 2; // C1
	 FF[3] = (1 + ksi / (1 - ksi)) * (1 + eta) / 2; // Q1

	 if (is_dyn) {
		 double pi = 3.14159265358979323;
		 double v_p = std::sqrt((Young * Poisson / ((1 + Poisson) * (1 - 2 * Poisson)) + Young / (2 * (1 + Poisson))) / density);
		 A = v_p / omega;
		 k = 2 * pi / A;

		 std::complex<double> i(0, 1);
		 double den = 1.0 - ksi;
		 std::complex<double> mult = sqrt(2.0 / den) * exp(k * A * i / 2.0) * exp(k * A * i / den);

		 std::transform(FF.begin(), FF.end(), FF.begin(), [&mult](std::complex<double> a) { return a * mult; });
	 }

	return FF;
}

double infQuad::gaussPoint(LocVar var, int i) {
	std::vector<std::vector<double>> gp = { { 0.3333333333, -1.0, -1.0, 0.3333333333 },
										  { -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926 },
										  { 0.0, 0.0, 0.0, 0.0 } };

	return gp[static_cast<int>(var)][i];
}

double infQuad::weight(LocVar var, int i) {
	std::vector<std::vector<double>> w = { { 1.5, 0.5, 0.5, 1.5 },
										 { 1.0, 1.0, 1.0, 1.0 } };

	return w[static_cast<int>(var)][i];
}

