#include "Quad.h"

std::vector<double> infQuad::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	 FF.resize(4);
	 FF[0] = (1 + ksi / (1 - ksi)) * (1 - eta) / 2; // Q0
	 FF[1] = (-ksi / (1 - ksi)) * (1 - eta) / 2; // C0
	 FF[2] = (-ksi / (1 - ksi)) * (1 + eta) / 2; // C1
	 FF[3] = (1 + ksi / (1 - ksi)) * (1 + eta) / 2; // Q1

	if (is_dyn) {
		double pi = 3.14159265358979323;
		double v_p = std::sqrt((Young * Poisson / ((1 + Poisson) * (1 - 2 * Poisson)) + Young / (2 * (1 + Poisson))) / density);
		A = v_p / (2 * pi * omega);
		double Young;
		k = 2 * pi / A;

		double numerator = 1;//sqrt(2.0 / (1.0 - ksi));
		double denominator = 1.0 - ksi;
		double mult = numerator * (cos(k * A / 2.0) * cos(k * A / denominator) - sin(k * A / 2.0) * sin(k * A / denominator));

		std::transform(FF.begin(), FF.end(), FF.begin(), [&mult](double a) { return a * mult; });
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
