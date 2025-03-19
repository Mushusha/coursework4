#include "QuadElements.h"

std::vector<double> infQuad::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF.resize(4);
	FF[0] = (1 + ksi / (1 - ksi)) * (1 - eta) / 2; // Q0
	FF[1] = (-ksi / (1 - ksi)) * (1 - eta) / 2; // C0
	FF[2] = (-ksi / (1 - ksi)) * (1 + eta) / 2; // C1
	FF[3] = (1 + ksi / (1 - ksi)) * (1 + eta) / 2; // Q1
	return FF;
}