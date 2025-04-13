#include "Tri.h"
#include "Data.h"

#include <math.h>

Eigen::MatrixXd Tri::C() {
	Eigen::Matrix3d C;
	for (int i = 0; i < 3; i++) {
		C(i, 0) = 1;
		C(i, 1) = x[i];
		C(i, 2) = y[i];
	}
	return C;
}

std::pair<int, int> Tri::edge_to_node(int edge) {
	if (edge != 2)
		return std::pair<int, int>(edge, edge + 1);
	else if (edge == 2)
		return std::pair<int, int>(edge, 0);
	else
		throw runtime_error("Error: wrong edge");
}

double Tri::len_edge(int edge) {
	std::pair<int, int> coords = edge_to_node(edge);
	return std::sqrt(std::pow((x[coords.first] - x[coords.second]), 2) + std::pow((y[coords.first] - y[coords.second]), 2));
}

Eigen::MatrixXd Tri::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXd B(3, 6);
	for (int i = 0; i < 3; i++) {
		B(0, 2 * i) = y[(1 + i) % 3] - y[(2 + i) % 3];
		B(0, 2 * i + 1) = 0;
		B(1, 2 * i) = 0;
		B(1, 2 * i + 1) = x[(2 + i) % 3] - x[(1 + i) % 3];
		B(2, 2 * i) = x[(2 + i) % 3] - x[(1 + i) % 3];
		B(2, 2 * i + 1) = y[(1 + i) % 3] - y[(2 + i) % 3];
	}
	B = B / std::abs(C().determinant());
	return B;
}

bool Tri::pointInElem(std::vector<double> point) {
	bool answer = true;
	for (int i = 0; i < 3; i++) {
		std::vector<double> a{ x[i], y[i] };
		std::vector<double> b{ x[(i + 1) % 3], y[(i + 1) % 3] };
		std::vector<double> c{ x[(i + 2) % 3], y[(i + 2) % 3] };
		if (a == point || b == point || c == point)
			return true;
		if (line(a, b, c) <= 0 && line(a, b, point) <= 0)
			answer &= true;
		else if (line(a, b, c) >= 0 && line(a, b, point) >= 0)
			answer &= true;
		else
			return false;
	}
	return answer;
}

Eigen::MatrixXd Tri::localK() {
	return B().transpose() * D * B() * std::abs(C().determinant() / 2);
}

std::vector<double> Tri::localF() {
	std::vector<double> F;
	F.resize(6);
	// l.first.first - edge, l.first.second - comp, l.second - value
	for (auto const& l: load) {
		std::pair<int, int> node = edge_to_node(l.first.first);
		F[2 * node.first + l.first.second] += l.second * len_edge(l.first.first) / 2;
		F[2 * node.second + l.first.second] += l.second * len_edge(l.first.first) / 2;
	}
	return F;
}

Eigen::MatrixXd Tri::localC() {
	Eigen::MatrixXd c(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			c(i, j) = (i == j) ? (std::abs(C().determinant() / 12)) : std::abs((C().determinant() / 24));
	return c;
}

std::vector<double> Tri::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(3);
	for (int i = 0; i < 3; i++)
		R[i] = value[i] * std::abs(C().determinant() / 6);
	return R;
}

Eigen::MatrixXd Tri::localM() {
	if (density == 0.0)
		throw runtime_error("Error: density is zero in element " + to_string(id));

	Eigen::MatrixXd m(6, 6);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			if ((i + j) % 2 == 1)
				m(i, j) = 0;
			else
				m(i, j) = (i == j) ? (std::abs(C().determinant() / 12)) : std::abs((C().determinant() / 24));
	return density * m;
}

std::vector<double> Tri::coordFF(double x0, double y0, double z0) {
	std::vector<double> coord = { x0, y0 };
	return coord;
}

double Tri::Volume() {
	return C().determinant();
}

std::vector<double> Tri::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF.resize(3);
	Eigen::Vector3d f = { 1, ksi, eta };
	Eigen::Vector3d ff = f.transpose() * C().inverse();
	FF = { ff(0), ff(1), ff(2) };
	//for (int i = 0; i < 3; i++) {
	//	Eigen::MatrixXd A(3, 3);
	//	A <<
	//		1, 1, 1,
	//		ksi, x[i], x[(i + 1) % 3],
	//		eta, y[i], y[(i + 1) % 3];
	//	FF[i] = std::abs(A.determinant() / C().determinant());
	//}
	return FF;
}

Eigen::MatrixXd Tri::gradFF(double ksi, double eta, double zeta) {
	return Eigen::MatrixXd();
}

Eigen::MatrixXd Tri::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXd j = Eigen::MatrixXd::Zero(2, 2);
	j(0, 0) = 1;
	j(1, 1) = 1;
	return j;
}

double Tri::gaussPoint(LocVar var, int i) {
	return 0.0;
}

double Tri::weight(LocVar var, int i) {
	return 0.0;
}

void Tri::set_pressure(int edge, double value) {
	std::pair<int, int> node = edge_to_node(edge);
	std::array<double, 2> comp;
	comp[0] = -y[node.first] + y[node.second];
	comp[1] = x[node.first] - x[node.second];

	if ((x[node.first] - x[node.second]) * (y[node.first] - y[3 - node.first - node.second]) -
		(y[node.first] - y[node.second]) * (x[node.first] - x[3 - node.first - node.second]) < 0)
		for (auto& i : comp)
			i *= -1;
	
	for (int i = 0; i < 2; i++) {
		std::pair <int, int> pair(edge, i);
		load.insert({ pair, -value * comp[i] / len_edge(edge)});
	}
}
