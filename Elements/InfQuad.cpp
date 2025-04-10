#include "Quad.h"

std::vector<double> infQuad::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF.resize(4);
	FF[0] = (1 + ksi / (1 - ksi)) * (1 - eta) / 2; // Q0
	FF[1] = (-ksi / (1 - ksi)) * (1 - eta) / 2; // C0
	FF[2] = (-ksi / (1 - ksi)) * (1 + eta) / 2; // C1
	FF[3] = (1 + ksi / (1 - ksi)) * (1 + eta) / 2; // Q1
	return FF;
}

std::vector<double> infQuadDyn::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF.resize(4);
	FF[0] = (1 + ksi / (1 - ksi)) * (1 - eta) / 2; // Q0
	FF[1] = (-ksi / (1 - ksi)) * (1 - eta) / 2; // C0
	FF[2] = (-ksi / (1 - ksi)) * (1 + eta) / 2; // C1
	FF[3] = (1 + ksi / (1 - ksi)) * (1 + eta) / 2; // Q1

	//std::complex<double> i(0, 1);
	//std::complex<double> numerator = sqrt(2.0 / (1.0 - ksi));
	//std::complex<double> denominator = 1.0 - ksi;
	//std::complex<double> mult = numerator * exp(i * k * A / 2.0) * exp(i * k * A / denominator);

	//std::transform(FF.begin(), FF.end(), FF.begin(), [&mult](std::complex<double> a) { return a * mult; });
	return FF;
}

double infQuadDyn::gaussPoint(LocVar var, int i) {
	return 0.0;
}
