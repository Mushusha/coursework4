#include "TriangleElements.h"
#include "Data.h"

#include <math.h>

Eigen::MatrixXd triElement::C() {
	Eigen::Matrix3d C;
	for (int i = 0; i < 3; i++) {
		C(i, 0) = 1;
		C(i, 1) = x[i];
		C(i, 2) = y[i];
	}
	return C;
}

std::pair<int, int> triElement::edge_to_node(int edge) {
	if (edge != 2)
		return std::pair<int, int>(edge, edge + 1);
	else if (edge == 2)
		return std::pair<int, int>(edge, 0);
	else
		std::cout << "Error: wrong edge" << endl;
}

double triElement::len_edge(int edge) {
	std::pair<int, int> coords = edge_to_node(edge);
	return std::sqrt(std::pow((x[coords.first] - x[coords.second]), 2) + std::pow((y[coords.first] - y[coords.second]), 2));
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
	B = B / std::abs(C().determinant());
	return B;
}

Eigen::MatrixXd triElement::localK() {
	return B().transpose() * planeStrainD() * B() * std::abs(C().determinant() / 2);
}

std::vector<double> triElement::localF() {
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

Eigen::MatrixXd triElement::localC() {
	Eigen::MatrixXd c(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			c(i, j) = (i == j) ? (C().determinant() / 12) : (C().determinant() / 24);
	return c;
}

std::vector<double> triElement::localR(double value) {
	std::vector<double> R;
	R.resize(3);
	for (int i = 0; i < 3; i++)
		R[i] = value * C().determinant() / 6;
	return R;
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

void triElement::set_pressure(int edge, double value) {
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
