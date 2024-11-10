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
	B = B / C().determinant() / 2;
	return B;
}

Eigen::MatrixXd triElement::locK() {
	return B().transpose() * twoMatrixD() * B() * C().determinant() / 2;
}

std::vector<double> triElement::locR() {
	std::vector<double> R;
	R.resize(6);
	// l.first.first - edge, l.first.second - comp, l.second - value
	for (auto const& l: load) {
		std::pair<int, int> nodes = edge_to_node(l.first.first);
		R[nodes.first + l.first.second] = l.second * len_edge(l.first.first) / 2;
		R[nodes.second + l.first.second] = l.second * len_edge(l.first.first) / 2;
	}
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