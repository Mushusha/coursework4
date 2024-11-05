#include "TriangleElements.h"
#include "Data.h"


Eigen::MatrixXd triElement::C() {
	Eigen::Matrix3d C;
	for (int i = 0; i < 3; i++) {
		C(i, 0) = 1;
		C(i, 1) = x[i];
		C(i, 2) = y[i];
	}
	return C;
}

Eigen::MatrixXd triElement::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXd B(3, 6);
	for (int i = 0; i < 3; i++) {
		B(0, 2 * i) = y[(1 + i) % 3] - y[(2 + i) % 3];
		B(0, 2 * i + 1) = 0;
		B(1, 2 * i) = 0;
		B(1, 2 * i + 1) = x[(2 + i) % 3] - x[(1 + i) % 3];
		B(2, 2 * i) = x[(2 + i) % 3] - x[(1 + i) % 3];
		B(2, 2 * i + 1) = y[(1 + i) % 3] - y[(2 + i) % 3];
	}
	B = B / C().determinant() / 2;
	return B;
}

Eigen::MatrixXd triElement::locK() {
	return B().transpose() * twoMatrixD() * B() * C().determinant() / 2;
}

std::vector<double> triElement::FF(double ksi, double eta, double zeta) {
	return std::vector<double>();
}

std::vector<std::vector<double>> triElement::gradFF(double ksi, double eta, double zeta) {
	return std::vector<std::vector<double>>();
}

Eigen::MatrixXd triElement::J(double ksi, double eta, double zeta) {
	return Eigen::MatrixXd();
}
