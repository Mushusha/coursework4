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

Eigen::MatrixXcd Tetra::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd B = Eigen::MatrixXd::Zero(6, 12);
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

Eigen::MatrixXcd Tetra::localK() {
	return B().transpose() * D * B() * C().determinant() / 6;
}

std::vector<double> Tetra::localF(double mult) {
	return std::vector<double>();
}

Eigen::MatrixXd Tetra::localC() {
	Eigen::MatrixXd c(4, 4);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			c(i, j) = (i == j) ? (std::abs(C().determinant() / 60)) : std::abs((C().determinant() / 120));
	return c;
}

std::vector<double> Tetra::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(4);
	for (int i = 0; i < 4; i++)
		R[i] = value[i] * std::abs(C().determinant() / 24);
	return R;
}

Eigen::MatrixXcd Tetra::localM() {
	if (density == 0.0)
		throw runtime_error("Error: density is zero in element " + to_string(id));

	Eigen::MatrixXd m(12, 12);
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
			if ((i + j) % 2 == 1)
				m(i, j) = 0;
			else
				m(i, j) = (i == j) ? (std::abs(C().determinant() / 60)) : std::abs((C().determinant() / 120));
	return density * m;
}

std::vector<double> Tetra::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}

double Tetra::Volume() {
	return C().determinant() / 6;
}

std::vector<std::complex<double>> Tetra::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
	FF.resize(4);
	Eigen::Vector4cd f = { 1, ksi, eta, zeta };
	Eigen::Vector4cd ff = f.transpose() * C().inverse();
	FF = { ff(0), ff(1), ff(2), ff(3)};
	return FF;
}

Eigen::MatrixXcd Tetra::gradFF(double ksi, double eta, double zeta) {
	return Eigen::MatrixXcd();
}

Eigen::MatrixXcd Tetra::J(double ksi, double eta, double zeta) {
	return Eigen::MatrixXcd();
}

double Tetra::gaussPoint(LocVar var, int i) {
	return 0.0;
}

double Tetra::weight(LocVar var, int i) {
	return 0.0;
}

void Tetra::set_pressure(int edge, double value) {
	//std::pair<int, int> node = edge_to_node(edge);
	//std::array<double, 2> comp;
	//comp[0] = -y[node.first] + y[node.second];
	//comp[1] = x[node.first] - x[node.second];

	//if ((x[node.first] - x[node.second]) * (y[node.first] - y[3 - node.first - node.second]) -
	//	(y[node.first] - y[node.second]) * (x[node.first] - x[3 - node.first - node.second]) < 0)
	//	for (auto& i : comp)
	//		i *= -1;

	//for (int i = 0; i < 2; i++) {
	//	std::pair <int, int> pair(edge, i);
	//	load.insert({ pair, -value * comp[i] / len_edge(edge) });
	//}
}
