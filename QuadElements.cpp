#include "QuadElements.h"

quadElement::quadElement(std::vector <int> _nodes, std::vector <float> _x, std::vector <float> _y) {
	for (int i = 0; i < _nodes.size(); i++) {
		nodes.push_back(_nodes[i]);
		x.push_back(_x[i]);
		y.push_back(_y[i]);
		z.push_back(0);
	}
	std::cout << "uhyu";
	locK = quadLocK(x, y, quadMapN);
}

quadInfElement::quadInfElement(std::vector <int> _nodes, std::vector <float> _x, std::vector <float> _y) {
	for (int i = 0; i < _nodes.size(); i++) {
		nodes.push_back(_nodes[i]);
		x.push_back(_x[i]);
		y.push_back(_y[i]);
		z.push_back(0);
	}
	locK = quadLocK(x, y, quadInfN);
}

float quadN1(float ksi, float eta) {
	return (1 - ksi) * (1 - eta) / 4;
}

float quadN2(float ksi, float eta) {
	return (1 + ksi) * (1 - eta) / 4;
}

float quadN3(float ksi, float eta) {
	return (1 + ksi) * (1 + eta) / 4;
}

float quadN4(float ksi, float eta) {
	return (1 - ksi) * (1 + eta) / 4;
}

float quadInfN(float ksi, float eta, std::vector <float> a) {
	return ((1 + eta) * (-ksi / (1 - ksi) * a[0] + (1 + ksi / (1 - ksi)) * a[1]) + (1 - eta) * (-ksi / (1 - ksi) * a[2] + (1 + ksi / (1 - ksi)) * a[3])) / 2;
}
// 0 - x_C, 1 - x_Q, 2 - x_C1, 3 - x_Q1

float quadMapN(float ksi, float eta, std::vector <float> a) {
	return quadN1(ksi, eta) * a[0] + quadN2(ksi, eta) * a[1] + quadN3(ksi, eta) * a[2] + quadN4(ksi, eta) * a[3];
}

//float shapeN(float ksi, float eta) {
//	return N1(ksi, eta) + N2(ksi, eta) + N3(ksi, eta) + N4(ksi, eta);
//}


Eigen::Matrix2f quadMatrixJ(float ksi, float eta, std::vector <float> x, std::vector <float> y, float (*N)(float, float, std::vector <float>)) {

	Eigen::Matrix2f J;
	float h = 0.01;
	J(0, 0) = (N(ksi + h, eta, x) - N(ksi - h, eta, x)) / (2 * h);
	J(0, 1) = (N(ksi, eta + h, x) - N(ksi, eta - h, x)) / (2 * h);
	J(1, 0) = (N(ksi + h, eta, y) - N(ksi - h, eta, y)) / (2 * h);
	J(1, 1) = (N(ksi, eta + h, y) - N(ksi, eta - h, y)) / (2 * h);
	return J;
}

float quadDksiN(float ksi, float eta, float (*N)(float, float)) {
	float h = 0.01;
	return (N(ksi + h, eta) - N(ksi - h, eta)) / (2 * h);
}

float quadDetaN(float ksi, float eta, float (*N)(float, float)) {
	float h = 0.01;
	return (N(ksi, eta + h) - N(ksi, eta - h)) / (2 * h);
}

Eigen::MatrixXf  quadMatrixB(float ksi, float eta, std::vector <float> x, std::vector <float> y, float (*funcN)(float, float, std::vector <float>)) {

	std::vector <arrFunc> N = { quadN1, quadN2, quadN3, quadN4 };
	Eigen::MatrixXf B = Eigen::MatrixXf::Zero(3, 8);
	Eigen::Matrix2f invJ;
	invJ = quadMatrixJ(ksi, eta, x, y, funcN).inverse();
	float h = 0.01;
	for (int i = 0; i < 4; i++) {
		B(0, 2 * i) = quadDksiN(ksi, eta, N[i]) * invJ(0, 0) + quadDetaN(ksi, eta, N[i]) * invJ(0, 1);
		B(1, 2 * i + 1) = quadDksiN(ksi, eta, N[i]) * invJ(1, 0) + quadDetaN(ksi, eta, N[i]) * invJ(1, 1);
		B(2, 2 * i) = quadDksiN(ksi, eta, N[i]) * invJ(1, 0) + quadDetaN(ksi, eta, N[i]) * invJ(1, 1);
		B(2, 2 * i + 1) = quadDksiN(ksi, eta, N[i]) * invJ(0, 0) + quadDetaN(ksi, eta, N[i]) * invJ(0, 1);
	}
	return B;
}

Eigen::MatrixXf quadLocK(std::vector <float> x, std::vector <float> y, float (*N)(float, float, std::vector <float>)) {

	Eigen::MatrixXf k = Eigen::MatrixXf::Zero(8, 8);
	std::vector <float> wieghtX = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <float> wieghtY = { 0.5774, 0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < wieghtX.size(); i++)
		k += quadMatrixB(wieghtX[i], wieghtY[i], x, y, N).transpose() * twoMatrixD() * quadMatrixB(wieghtX[i], wieghtY[i], x, y, N) * abs(quadMatrixJ(wieghtX[i], wieghtY[i], x, y, N).determinant());
	return k;
}