#include "Tetra.h"


Eigen::MatrixXd Tetra::C() {

	Eigen::Matrix4d C;
	for (int i = 0; i < 4; i++) {
		C(i, 0) = 1;
		C(i, 1) = x[i];
		C(i, 2) = y[i];
		C(i, 3) = z[i];

	}
	return C;
}

Eigen::MatrixXd Tetra::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 12);
	Eigen::Matrix4d coef;
	Eigen::Matrix3d a, b, c, d;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			a(j, 0) = x[(i + j + 1) % 4];
			a(j, 1) = y[(i + j + 1) % 4];
			a(j, 2) = z[(i + j + 1) % 4];
	
			b(j, 0) = 1;
			b(j, 1) = y[(i + j + 1) % 4];
			b(j, 2) = z[(i + j + 1) % 4];
	
			c(j, 0) = x[(i + j + 1) % 4];
			c(j, 1) = 1;
			c(j, 2) = z[(i + j + 1) % 4];
	
			d(j, 0) = x[(i + j + 1) % 4];
			d(j, 1) = y[(i + j + 1) % 4];
			d(j, 2) = 1;
		}
		coef(0, i) = a.determinant();
		coef(1, i) = -1 * b.determinant();
		coef(2, i) = -1 * c.determinant();
		coef(3, i) = -1 * d.determinant();
	}

	for (int i = 0; i < 4; i++) {
		B(0, i * 3) = coef(i, 1);
		B(1, i * 3 + 1) = coef(i, 2);
		B(2, i * 3 + 2) = coef(i, 3);
		B(3, i * 3) = coef(i, 2);
		B(3, i * 3 + 1) = coef(i, 1);
		B(4, i * 3 + 1) = coef(i, 3);
		B(4, i * 3 + 2) = coef(i, 2);
		B(5, i * 3) = coef(i, 3);
		B(5, i * 3 + 2) = coef(i, 1);
	}

	return B / (6 * C().determinant());
}

bool Tetra::pointInElem(std::vector<double> point) {
	return false;
}

Eigen::MatrixXd Tetra::localK() {
	return B().transpose() * D * B() * C().determinant() / 6;
}

std::vector<double> Tetra::localF() {
	return std::vector<double>();
}

Eigen::MatrixXd Tetra::localC() {
	return Eigen::MatrixXd();
}

std::vector<double> Tetra::localR(std::vector<double> value) {
	return std::vector<double>();
}

Eigen::MatrixXd Tetra::localM() {
	return Eigen::MatrixXd();
}

std::vector<double> Tetra::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}

double Tetra::Volume() {
	return 0.0;
}

std::vector<double> Tetra::FF(double ksi, double eta, double zeta) {
	return std::vector<double>();
}

Eigen::MatrixXd Tetra::gradFF(double ksi, double eta, double zeta) {
	return Eigen::MatrixXd();
}

Eigen::MatrixXd Tetra::J(double ksi, double eta, double zeta) {
	return Eigen::MatrixXd();
}

double Tetra::gaussPoint(LocVar var, int i) {
	return 0.0;
}

double Tetra::weight(LocVar var, int i) {
	return 0.0;
}

void Tetra::set_pressure(int edge, double value) {
}
