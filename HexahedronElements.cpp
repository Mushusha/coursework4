#include "HexahedronElements.h"

float hexaN1(float ksi, float eta, float zeta) {
	return (1 - ksi) * (1 - eta) * (1 - zeta) / 8;
}

float hexaN2(float ksi, float eta, float zeta) {
	return (1 + ksi) * (1 - eta) * (1 - zeta) / 8;
}

float hexaN3(float ksi, float eta, float zeta) {
	return (1 + ksi) * (1 + eta) * (1 - zeta) / 8;
}

float hexaN4(float ksi, float eta, float zeta) {
	return (1 - ksi) * (1 + eta) * (1 - zeta) / 8;
}

float hexaN5(float ksi, float eta, float zeta) {
	return (1 - ksi) * (1 - eta) * (1 + zeta) / 8;
}

float hexaN6(float ksi, float eta, float zeta) {
	return (1 + ksi) * (1 - eta) * (1 + zeta) / 8;
}

float hexaN7(float ksi, float eta, float zeta) {
	return (1 + ksi) * (1 + eta) * (1 + zeta) / 8;
}

float hexaN8(float ksi, float eta, float zeta) {
	return (1 - ksi) * (1 + eta) * (1 + zeta) / 8;
}

float hexaInfN(float ksi, float eta, float zeta, std::vector <float> a) {
	return ((-ksi / (1 - ksi)) * ((1 - eta) * (1 - zeta) * a[0] + (1 + eta) * (1 - zeta) * a[1] +
		(1 + eta) * (1 + zeta) * a[2] + (1 - eta) * (1 + zeta) * a[3])
		+ (1 + ksi / (1 - ksi)) * ((1 - eta) * (1 - zeta) * a[4] + (1 + eta) * (1 - zeta) * a[5] +
			(1 + eta) * (1 + zeta) * a[6] + (1 - eta) * (1 + zeta) * a[7])) / 4;
}
// 0 - x_C, 1 - x_C1, 2 - x_C2, 3 - x_C3, 4 - x_Q, 5 - x_Q1, 6 = x_Q2, 7 - x_Q3

float hexaMapN(float ksi, float eta, float zeta, std::vector <float> a) {
	return hexaN1(ksi, eta, zeta) * a[0] + hexaN2(ksi, eta, zeta) * a[1] + hexaN3(ksi, eta, zeta) * a[2] + hexaN4(ksi, eta, zeta) * a[3] +
		hexaN5(ksi, eta, zeta) * a[4] + hexaN6(ksi, eta, zeta) * a[5] + hexaN7(ksi, eta, zeta) * a[6] + hexaN8(ksi, eta, zeta) * a[7];
}

//float shapeN(float ksi, float eta) {
//	return N1(ksi, eta) + N2(ksi, eta) + N3(ksi, eta) + N4(ksi, eta);
//}

Eigen::Matrix3f hexaMatrixJ(float ksi, float eta, float zeta, std::vector <float> x, std::vector <float> y, std::vector <float> z, float (*N)(float, float, float, std::vector <float>)) {

	Eigen::Matrix3f J;
	float h = 0.01;

	J(0, 0) = (N(ksi + h, eta, zeta, x) - N(ksi - h, eta, zeta, x)) / (2 * h);
	J(0, 1) = (N(ksi, eta + h, zeta, x) - N(ksi, eta - h, zeta, x)) / (2 * h);
	J(0, 2) = (N(ksi, eta, zeta + h, x) - N(ksi, eta, zeta - h, x)) / (2 * h);
	J(1, 0) = (N(ksi + h, eta, zeta, y) - N(ksi - h, eta, zeta, y)) / (2 * h);
	J(1, 1) = (N(ksi, eta + h, zeta, y) - N(ksi, eta - h, zeta, y)) / (2 * h);
	J(1, 2) = (N(ksi, eta, zeta + h, y) - N(ksi, eta, zeta - h, y)) / (2 * h);
	J(2, 0) = (N(ksi + h, eta, zeta, z) - N(ksi - h, eta, zeta, z)) / (2 * h);
	J(2, 1) = (N(ksi, eta + h, zeta, z) - N(ksi, eta - h, zeta, z)) / (2 * h);
	J(2, 2) = (N(ksi, eta, zeta + h, z) - N(ksi, eta, zeta - h, z)) / (2 * h);

	return J;
}

float hexaDksiN(float ksi, float eta, float zeta, float (*N)(float, float, float)) {
	float h = 0.01;
	return (N(ksi + h, eta, zeta) - N(ksi - h, eta, zeta)) / (2 * h);
}

float hexaDetaN(float ksi, float eta, float zeta, float (*N)(float, float, float)) {
	float h = 0.01;
	return (N(ksi, eta + h, zeta) - N(ksi, eta - h, zeta)) / (2 * h);
}

float hexaDzetaN(float ksi, float eta, float zeta, float (*N)(float, float, float)) {
	float h = 0.01;
	return (N(ksi, eta, zeta + h) - N(ksi, eta, zeta - h)) / (2 * h);
}

Eigen::MatrixXf  hexaMatrixB(float ksi, float eta, float zeta, std::vector <float> x, std::vector <float> y, std::vector <float> z, float (*funcN)(float, float, float, std::vector <float>)) {

	std::vector <hexaArrFunc> N = { hexaN1, hexaN2, hexaN3, hexaN4, hexaN5, hexaN6, hexaN7, hexaN8 };
	Eigen::MatrixXf B = Eigen::MatrixXf::Zero(6, 24);
	Eigen::Matrix3f invJ;
	invJ = hexaMatrixJ(ksi, eta, zeta, x, y, z, funcN).inverse();
	float h = 0.01;
	for (int i = 0; i < 8; i++) {
		B(0, 3 * i) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(0, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(0, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(0, 2);
		B(1, 3 * i + 1) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(1, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(1, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(1, 2);
		B(2, 3 * i + 2) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(2, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(2, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(2, 2);
		B(3, 3 * i) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(1, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(1, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(1, 2);
		B(3, 3 * i + 1) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(0, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(0, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(0, 2);
		B(4, 3 * i + 1) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(2, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(2, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(2, 2);
		B(4, 3 * i + 2) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(1, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(1, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(1, 2);
		B(5, 3 * i) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(2, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(2, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(2, 2);
		B(5, 3 * i + 2) = hexaDksiN(ksi, eta, zeta, N[i]) * invJ(0, 0) + hexaDetaN(ksi, eta, zeta, N[i]) * invJ(0, 1) + hexaDzetaN(ksi, eta, zeta, N[i]) * invJ(0, 2);
	}
	return B;
}

Eigen::MatrixXf hexaLocK(std::vector <float> x, std::vector <float> y, std::vector <float> z, float (*N)(float, float, float, std::vector <float>)) {

	Eigen::MatrixXf k = Eigen::MatrixXf::Zero(24, 24);
	std::vector <float> wieghtX = { 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <float> wieghtY = { 0.5774, 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774 };
	std::vector <float> wieghtZ = { 0.5774, 0.5774, 0.5774, 0.5774, -0.5774, -0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < wieghtX.size(); i++)
		k += hexaMatrixB(wieghtX[i], wieghtY[i], wieghtZ[i], x, y, z, N).transpose() * threeMatrixD() * hexaMatrixB(wieghtX[i], wieghtY[i], wieghtZ[i], x, y, z, N) * abs(hexaMatrixJ(wieghtX[i], wieghtY[i], wieghtZ[i], x, y, z, N).determinant());
	return k;
}
