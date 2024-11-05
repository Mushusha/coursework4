#include "QuadElements.h"

std::vector<double> quadElement::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF[0] = (1 - ksi) * (1 - eta) / 4;
	FF[1] = (1 + ksi) * (1 - eta) / 4;
	FF[2] = (1 + ksi) * (1 + eta) / 4;
	FF[3] = (1 - ksi) * (1 + eta) / 4;
	return FF;
}

std::vector<std::vector<double>> quadElement::gradFF(double ksi, double eta, double zeta) {
	std::vector <std::vector <double>> gradFF;
	double h = 0.01;
	for (int i = 0; i < 4; i++) {
		gradFF[LocVar::KSI][i] = (FF(ksi + h, eta)[i] - FF(ksi - h, eta)[i]) / (2 * h);
		gradFF[LocVar::ETA][i] = (FF(ksi, eta + h)[i] - FF(ksi, eta - h)[i]) / (2 * h);
	}
	return gradFF;
}

Eigen::MatrixXd quadElement::J(double ksi, double eta, double zeta) {
	Eigen::Matrix2d J;
	double h = 0.01;
	for (int i = 0; i < 4; i++) {
		J(0, 0) += gradFF(ksi, eta)[LocVar::KSI][i] * x[i];
		J(0, 1) += gradFF(ksi, eta)[LocVar::ETA][i] * x[i];
		J(1, 0) += gradFF(ksi, eta)[LocVar::KSI][i] * y[i];
		J(1, 1) += gradFF(ksi, eta)[LocVar::ETA][i] * y[i];
	}
	return J;
}

Eigen::MatrixXd  quadElement::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
	Eigen::Matrix2d invJ;
	invJ = J(ksi, eta).inverse();
	double h = 0.01;
	for (int i = 0; i < 4; i++) {
		B(0, 2 * i) = gradFF(ksi, eta)[LocVar::KSI][i] * invJ(0, 0) + gradFF(ksi, eta)[LocVar::ETA][i] * invJ(0, 1);
		B(1, 2 * i + 1) = gradFF(ksi, eta)[LocVar::KSI][i] * invJ(1, 0) + gradFF(ksi, eta)[LocVar::ETA][i] * invJ(1, 1);
		B(2, 2 * i) = gradFF(ksi, eta)[LocVar::KSI][i] * invJ(1, 0) + gradFF(ksi, eta)[LocVar::ETA][i] * invJ(1, 1);
		B(2, 2 * i + 1) = gradFF(ksi, eta)[LocVar::KSI][i] * invJ(0, 0) + gradFF(ksi, eta)[LocVar::ETA][i] * invJ(0, 1);
	}
	return B;
}

Eigen::MatrixXd quadElement::locK() {

	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(8, 8);
	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < 4; i++)
		k += B(ksi[i], eta[i]).transpose() * twoMatrixD() * B(ksi[i], eta[i]) * std::abs(J(eta[i], ksi[i]).determinant());
	return k;
}

//double quadInfN(double ksi, double eta, std::vector <double> a) {
//	return ((1 + eta) * (-ksi / (1 - ksi) * a[0] + (1 + ksi / (1 - ksi)) * a[1]) + (1 - eta) * (-ksi / (1 - ksi) * a[2] + (1 + ksi / (1 - ksi)) * a[3])) / 2;
//}
// 0 - x_C, 1 - x_Q, 2 - x_C1, 3 - x_Q1
