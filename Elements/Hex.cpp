#include "Hex.h"


std::vector<std::complex<double>> Hex::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
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

Eigen::MatrixXcd Hex::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd gradFF;
	//for (int i = 0; i < 8; i++) {
	//	gradFF[KSI][i] = (FF(ksi + h, eta, zeta)[i] - FF(ksi - h, eta, zeta)[i]) / (2 * h);
	//	gradFF[ETA][i] = (FF(ksi, eta + h, zeta)[i] - FF(ksi, eta - h, zeta)[i]) / (2 * h);
	//	gradFF[ZETA][i] = (FF(ksi, eta, zeta + h)[i] - FF(ksi, eta, zeta - h)[i]) / (2 * h);
	//}
	return gradFF;
}

Eigen::MatrixXcd Hex::J(double ksi, double eta, double zeta) {
	Eigen::Matrix3cd J;
	double h = 0.01;
	for (int i = 0; i < 8; i++) {
		//J(0, 0) = gradFF(ksi, eta, zeta)[KSI][i] * x[i];
		//J(0, 1) = gradFF(ksi, eta, zeta)[ETA][i] * x[i];
		//J(0, 2) = gradFF(ksi, eta, zeta)[ZETA][i] * x[i];
		//J(1, 0) = gradFF(ksi, eta, zeta)[KSI][i] * y[i];
		//J(1, 1) = gradFF(ksi, eta, zeta)[ETA][i] * y[i];
		//J(1, 2) = gradFF(ksi, eta, zeta)[ZETA][i] * y[i];
		//J(2, 0) = gradFF(ksi, eta, zeta)[KSI][i] * z[i];
		//J(2, 1) = gradFF(ksi, eta, zeta)[ETA][i] * z[i];
		//J(2, 2) = gradFF(ksi, eta, zeta)[ZETA][i] * z[i];
	}
	return J;
}

double Hex::gaussPoint(LocVar var, int i) {
	//std::vector <double> ksi = { -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926 };
	//std::vector <double> eta = { -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926 };
	//std::vector <double> zeta = { -0.57735026918926, -0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, 0.57735026918926, 0.57735026918926 };

	return 0.0;
}

double Hex::weight(LocVar var, int i) {
	return 0.0;
}

Eigen::MatrixXcd Hex::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(6, 24);
	Eigen::Matrix3cd invJ;
	invJ = J(ksi, eta, zeta).inverse();
	double h = 0.01;
	for (int i = 0; i < 8; i++) {
		//B(0, 3 * i) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(0, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(0, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(0, 2);
		//B(1, 3 * i + 1) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(1, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(1, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(1, 2);
		//B(2, 3 * i + 2) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(2, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(2, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(2, 2);
		//B(3, 3 * i) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(1, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(1, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(1, 2);
		//B(3, 3 * i + 1) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(0, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(0, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(0, 2);
		//B(4, 3 * i + 1) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(2, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(2, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(2, 2);
		//B(4, 3 * i + 2) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(1, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(1, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(1, 2);
		//B(5, 3 * i) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(2, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(2, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(2, 2);
		//B(5, 3 * i + 2) = gradFF(ksi, eta, zeta)[KSI][i] * invJ(0, 0) + gradFF(ksi, eta, zeta)[ETA][i] * invJ(0, 1) + gradFF(ksi, eta, zeta)[ZETA][i] * invJ(0, 2);
	}
	return B;
}

bool Hex::pointInElem(std::vector<double> point) {
	return false;
}

void Hex::set_pressure(int edge, double value) {
}

Eigen::MatrixXcd Hex::localK() {
	Eigen::MatrixXcd k = Eigen::MatrixXd::Zero(24, 24);
	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774 };
	std::vector <double> zeta = { 0.5774, 0.5774, 0.5774, 0.5774, -0.5774, -0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < 8; i++) // fix
		k += B(ksi[i], eta[i], zeta[i]).transpose() * D * B(ksi[i], eta[i], zeta[i]) * std::abs(J(ksi[i], eta[i], zeta[i]).determinant());
	return k;
}

std::vector<double> Hex::localF(double mult) {
	return std::vector<double>();
}

Eigen::MatrixXd Hex::localC() {
	return Eigen::MatrixXd();
}

std::vector<double> Hex::localR(std::vector<double> value) {
	return std::vector<double>();
}

Eigen::MatrixXcd Hex::localM() {
	return Eigen::MatrixXcd();
}

std::vector<double> Hex::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}

double Hex::Volume() {
	return 0.0;
}

//double hexInfN(double ksi, double eta, double zeta, std::vector <double> a) {
//	return ((-ksi / (1 - ksi)) * ((1 - eta) * (1 - zeta) * a[0] + (1 + eta) * (1 - zeta) * a[1] +
//		(1 + eta) * (1 + zeta) * a[2] + (1 - eta) * (1 + zeta) * a[3])
//		+ (1 + ksi / (1 - ksi)) * ((1 - eta) * (1 - zeta) * a[4] + (1 + eta) * (1 - zeta) * a[5] +
//			(1 + eta) * (1 + zeta) * a[6] + (1 - eta) * (1 + zeta) * a[7])) / 4;
//}
// 0 - x_C, 1 - x_C1, 2 - x_C2, 3 - x_C3, 4 - x_Q, 5 - x_Q1, 6 = x_Q2, 7 - x_Q3
