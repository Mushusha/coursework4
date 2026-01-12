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
	Eigen::Matrix3d b, c, d;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			int ind = (i + j + 1) % 4;

			b(j, 0) = 1;
			b(j, 1) = y[ind];
			b(j, 2) = z[ind];

			c(j, 0) = x[ind];
			c(j, 1) = 1;
			c(j, 2) = z[ind];

			d(j, 0) = x[ind];
			d(j, 1) = y[ind];
			d(j, 2) = 1;
		}
		int mult = (i % 2 == 0) ? -1 : 1;
		coef(1, i) = mult * b.determinant();
		coef(2, i) = mult * c.determinant();
		coef(3, i) = mult * d.determinant();
	}

	for (int i = 0; i < 4; i++) {
		B(0, i * 3) = coef(1, i);
		B(1, i * 3 + 1) = coef(2, i);
		B(2, i * 3 + 2) = coef(3, i);
		B(3, i * 3) = coef(2, i);
		B(3, i * 3 + 1) = coef(1, i);
		B(4, i * 3 + 1) = coef(3, i);
		B(4, i * 3 + 2) = coef(2, i);
		B(5, i * 3) = coef(3, i);
		B(5, i * 3 + 2) = coef(1, i);
	}

	return B / std::abs(C().determinant());
}

Eigen::MatrixXcd Tetra::localK() {
	return B().transpose() * D * B() * std::abs(C().determinant()) / 6;
}

std::vector<double> Tetra::localF(double mult) {
	std::vector<double> F;
	F.resize(12);

	// l.first.first - edge, l.first.second - comp, l.second - value
	for (auto const& l : load) {
		std::vector<int> node = edge_to_node(l.first.first);
		for (int i = 0; i < 3; i++)
			F[3 * node[i] + l.first.second] += mult * l.second / 3;
	}
	return F;
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

std::vector<std::complex<double>> Tetra::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
	FF.resize(4);
	Eigen::Vector4cd f = { 1, ksi, eta, zeta };
	Eigen::Vector4cd ff = f.transpose() * C().inverse();
	FF = { ff(0), ff(1), ff(2), ff(3)};
	return FF;
}

double Tetra::Volume() {
	return C().determinant() / 6;
}

std::vector<int> Tetra::edge_to_node(int edge) {
	switch (edge) {
	case 0:
		return { 0, 2, 3 };
	case 1:
		return { 2, 3, 1 };
	case 2:
		return { 0, 1, 3 };
	case 3:
		return { 1, 0, 2 };
	default:
		throw runtime_error("Error: wrong edge");
	}
}

std::array<double, 3> Tetra::normal(int edge) {
	std::vector<int> node = edge_to_node(edge);

	double v1x = x[node[1]] - x[node[0]];
	double v1y = y[node[1]] - y[node[0]];
	double v1z = z[node[1]] - z[node[0]];

	double v2x = x[node[2]] - x[node[0]];
	double v2y = y[node[2]] - y[node[0]];
	double v2z = z[node[2]] - z[node[0]];

	double nx = v1y * v2z - v1z * v2y;
	double ny = v1z * v2x - v1x * v2z;
	double nz = v1x * v2y - v1y * v2x;

	return { nx, ny, nz };
}

double Tetra::area_edge(int edge) {
	std::array<double, 3> n = normal(edge);
	return std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / 2;
}

void Tetra::set_pressure(int edge, double value) {
	std::vector<int> node = edge_to_node(edge);
	std::array<double, 3> comp;

	double area = area_edge(edge);

	std::array<double, 3> n = normal(edge);
	double len_norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

	for (auto& i : n)
		i /= len_norm;
	
	for (int i = 0; i < 3; i++) {
		std::pair <int, int> pair(edge, i);
		load.insert({ pair, -value * n[i] * area });
	}
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

bool Tetra::pointInElem(std::vector<double> point) {
	return false;
}

std::vector<double> Tetra::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}

