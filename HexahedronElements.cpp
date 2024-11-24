#include "HexahedronElements.h"


std::vector<double> hexElement::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF[0] = (1 - ksi) * (1 - eta) * (1 - zeta) / 8;
	FF[1] = (1 + ksi) * (1 - eta) * (1 - zeta) / 8;
	FF[2] = (1 + ksi) * (1 + eta) * (1 - zeta) / 8;
	FF[3] = (1 - ksi) * (1 + eta) * (1 - zeta) / 8;
	FF[4] = (1 - ksi) * (1 - eta) * (1 + zeta) / 8;
	FF[5] = (1 + ksi) * (1 - eta) * (1 + zeta) / 8;
	FF[6] = (1 + ksi) * (1 + eta) * (1 + zeta) / 8;
	FF[7] = (1 - ksi) * (1 + eta) * (1 + zeta) / 8;
	return FF;
}

std::vector <std::vector <double>> hexElement::gradFF(double ksi, double eta, double zeta) {
	std::vector <std::vector <double>> gradFF;
	double h = 0.01;
	for (int i = 0; i < 8; i++) {
		gradFF[LocVar::KSI][i] = (FF(ksi + h, eta, zeta)[i] - FF(ksi - h, eta, zeta)[i]) / (2 * h);
		gradFF[LocVar::ETA][i] = (FF(ksi, eta + h, zeta)[i] - FF(ksi, eta - h, zeta)[i]) / (2 * h);
		gradFF[LocVar::ZETA][i] = (FF(ksi, eta, zeta + h)[i] - FF(ksi, eta, zeta - h)[i]) / (2 * h);
	}
	return gradFF;
}

Eigen::MatrixXd hexElement::J(double ksi, double eta, double zeta) {
	Eigen::Matrix3d J;
	double h = 0.01;
	for (int i = 0; i < 8; i++) {
		J(0, 0) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * x[i];
		J(0, 1) = gradFF(ksi, eta, zeta)[LocVar::ETA][i] * x[i];
		J(0, 2) = gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * x[i];
		J(1, 0) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * y[i];
		J(1, 1) = gradFF(ksi, eta, zeta)[LocVar::ETA][i] * y[i];
		J(1, 2) = gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * y[i];
		J(2, 0) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * z[i];
		J(2, 1) = gradFF(ksi, eta, zeta)[LocVar::ETA][i] * z[i];
		J(2, 2) = gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * z[i];
	}
	return J;
}

Eigen::MatrixXd hexElement::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 24);
	Eigen::Matrix3d invJ;
	invJ = J(ksi, eta, zeta).inverse();
	double h = 0.01;
	for (int i = 0; i < 8; i++) {
		B(0, 3 * i) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(0, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(0, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(0, 2);
		B(1, 3 * i + 1) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(1, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(1, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(1, 2);
		B(2, 3 * i + 2) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(2, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(2, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(2, 2);
		B(3, 3 * i) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(1, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(1, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(1, 2);
		B(3, 3 * i + 1) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(0, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(0, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(0, 2);
		B(4, 3 * i + 1) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(2, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(2, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(2, 2);
		B(4, 3 * i + 2) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(1, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(1, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(1, 2);
		B(5, 3 * i) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(2, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(2, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(2, 2);
		B(5, 3 * i + 2) = gradFF(ksi, eta, zeta)[LocVar::KSI][i] * invJ(0, 0) + gradFF(ksi, eta, zeta)[LocVar::ETA][i] * invJ(0, 1) + gradFF(ksi, eta, zeta)[LocVar::ZETA][i] * invJ(0, 2);
	}
	return B;
}

void hexElement::set_pressure(int edge, double value) {
}

Eigen::MatrixXd hexElement::localK() {
	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(24, 24);
	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774 };
	std::vector <double> zeta = { 0.5774, 0.5774, 0.5774, 0.5774, -0.5774, -0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < 8; i++)
		k += B(ksi[i], eta[i], zeta[i]).transpose() * threeMatrixD() * B(ksi[i], eta[i], zeta[i]) * std::abs(J(ksi[i], eta[i], zeta[i]).determinant());
	return k;
}

std::vector<double> hexElement::localF() {
	return std::vector<double>();
}

Eigen::MatrixXd hexElement::localC() {
	return Eigen::MatrixXd();
}

std::vector<double> hexElement::localR(double value) {
	return std::vector<double>();
}

//double hexInfN(double ksi, double eta, double zeta, std::vector <double> a) {
//	return ((-ksi / (1 - ksi)) * ((1 - eta) * (1 - zeta) * a[0] + (1 + eta) * (1 - zeta) * a[1] +
//		(1 + eta) * (1 + zeta) * a[2] + (1 - eta) * (1 + zeta) * a[3])
//		+ (1 + ksi / (1 - ksi)) * ((1 - eta) * (1 - zeta) * a[4] + (1 + eta) * (1 - zeta) * a[5] +
//			(1 + eta) * (1 + zeta) * a[6] + (1 - eta) * (1 + zeta) * a[7])) / 4;
//}
// 0 - x_C, 1 - x_C1, 2 - x_C2, 3 - x_C3, 4 - x_Q, 5 - x_Q1, 6 = x_Q2, 7 - x_Q3
