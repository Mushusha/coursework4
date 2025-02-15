#include "QuadElements.h"

std::vector<double> infQuadElement::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF.resize(4);
	FF[0] = (1 + ksi / (1 - ksi)) * (1 - eta) / 2;
	FF[1] = (-ksi / (1 - ksi)) * (1 - eta) / 2;
	FF[2] = (-ksi / (1 - ksi)) * (1 + eta) / 2;
	FF[3] = (1 + ksi / (1 - ksi)) * (1 + eta) / 2;
	return FF;
}

//double quadInfN(double ksi, double eta, std::vector <double> a) {
//	return ((1 + eta) * (-ksi / (1 - ksi) * a[0] + (1 + ksi / (1 - ksi)) * a[1]) + (1 - eta) * (-ksi / (1 - ksi) * a[2] + (1 + ksi / (1 - ksi)) * a[3])) / 2;
//}
// 0 - x_C, 1 - x_Q, 2 - x_C1, 3 - x_Q1
