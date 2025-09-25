#include "Wedge.h"


std::vector<std::complex<double>> Wedge::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
	FF.resize(6);
	double L1 = (1 - ksi - eta) / 2;
	double L2 = ksi; // (1 + ksi) / 2;
	double L3 = eta; // (1 + eta) / 2;

	FF[0] = L1 * (1 - zeta) / 2;
	FF[1] = L2 * (1 - zeta) / 2;
	FF[2] = L3 * (1 - zeta) / 2;
	FF[3] = L1 * (1 + zeta) / 2;
	FF[4] = L2 * (1 + zeta) / 2;
	FF[5] = L3 * (1 + zeta) / 2;
	return FF;
}

Eigen::MatrixXcd Wedge::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd gradFF = Eigen::MatrixXcd::Zero(3, 6);
	double h = 0.01;
	for (int i = 0; i < 6; i++) {
		gradFF(KSI, i) = (FF(ksi + h, eta, zeta)[i] - FF(ksi - h, eta, zeta)[i]) / (2 * h);
		gradFF(ETA, i) = (FF(ksi, eta + h, zeta)[i] - FF(ksi, eta - h, zeta)[i]) / (2 * h);
		gradFF(ZETA, i) = (FF(ksi, eta, zeta + h)[i] - FF(ksi, eta, zeta - h)[i]) / (2 * h);
	}
	return gradFF;
}

Eigen::MatrixXcd Wedge::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd J = Eigen::MatrixXcd::Zero(3, 3);

	for (int i = 0; i < 6; i++) {
		J(0, 0) += gradFF(ksi, eta, zeta)(KSI, i) * x[i];
		J(0, 1) += gradFF(ksi, eta, zeta)(ETA, i) * x[i];
		J(0, 2) += gradFF(ksi, eta, zeta)(ZETA, i) * x[i];
		J(1, 0) += gradFF(ksi, eta, zeta)(KSI, i) * y[i];
		J(1, 1) += gradFF(ksi, eta, zeta)(ETA, i) * y[i];
		J(1, 2) += gradFF(ksi, eta, zeta)(ZETA, i) * y[i];
		J(2, 0) += gradFF(ksi, eta, zeta)(KSI, i) * z[i];
		J(2, 1) += gradFF(ksi, eta, zeta)(ETA, i) * z[i];
		J(2, 2) += gradFF(ksi, eta, zeta)(ZETA, i) * z[i];
	}
	return J;
}

Eigen::MatrixXcd Wedge::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(6, 18);
	Eigen::Matrix3cd invJ;
	invJ = J(ksi, eta, zeta).inverse();
	Eigen::MatrixXcd dN = invJ.transpose() * gradFF(ksi, eta, zeta);

	for (int i = 0; i < 6; i++) {
		B(0, 3 * i) = dN(X, i);
		B(1, 3 * i + 1) = dN(Y, i);
		B(2, 3 * i + 2) = dN(Z, i);
		B(3, 3 * i) = dN(Y, i);
		B(3, 3 * i + 1) = dN(X, i);
		B(4, 3 * i + 1) = dN(Z, i);
		B(4, 3 * i + 2) = dN(Y, i);
		B(5, 3 * i) = dN(Z, i);
		B(5, 3 * i + 2) = dN(X, i);
	}
	return B;
}

Eigen::MatrixXcd Wedge::localK() {
	Eigen::MatrixXcd k = Eigen::MatrixXcd::Zero(18, 18);

	for (int gp = 0; gp < 6; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
		double zeta = gaussPoint(ZETA, gp);

		k += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) * B(ksi, eta, zeta).transpose() *
			D * B(ksi, eta, zeta) * std::abs(J(ksi, eta, zeta).determinant());
	}
	return k;
}

std::vector<double> Wedge::localF(double mult) {
	std::vector<double> F;
	F.resize(18);

	// l.first.first - edge, l.first.second - comp, l.second - value
	if (load.size() != 0)
		for (auto const& l : load) {
			std::vector<int> node = edge_to_node(l.first.first);
			for (int i = 0; i < node.size(); i++)
				F[3 * node[i] + l.first.second] += mult * l.second / 4;
		}
	return F;
}

Eigen::MatrixXd Wedge::localC() {
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(6, 6);

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			for (int gp = 0; gp < 6; gp++) {
				double ksi = gaussPoint(KSI, gp);
				double eta = gaussPoint(ETA, gp);
				double zeta = gaussPoint(ZETA, gp);

				c(i, j) += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) * FF(ksi, eta, zeta)[i].real() *
					FF(ksi, eta, zeta)[j].real() * std::abs(J(ksi, eta, zeta).determinant());
			}
	return c;
}

std::vector<double> Wedge::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(6);

	for (int i = 0; i < 6; i++)
		for (int gp = 0; gp < 6; gp++) {
			double ksi = gaussPoint(KSI, gp);
			double eta = gaussPoint(ETA, gp);
			double zeta = gaussPoint(ZETA, gp);

			R[i] += value[i] * weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
				FF(ksi, eta, zeta)[i].real() * std::abs(J(ksi, eta, zeta).determinant());
		}
	return R;
}

Eigen::MatrixXcd Wedge::localM() {
	if (density == 0.0)
		throw runtime_error("Error: density is zero in element " + to_string(id));

	Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(18, 18);

	for (size_t gp = 0; gp < 6; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
		double zeta = gaussPoint(ZETA, gp);

		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 6; j++) {
				complex<double> M_ij = weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
					FF(ksi, eta, zeta)[i] * FF(ksi, eta, zeta)[j] * std::abs(J(ksi, eta, zeta).determinant());

				m(3 * i, 3 * j) += M_ij;
				m(3 * i + 1, 3 * j + 1) += M_ij;
				m(3 * i + 2, 3 * j + 2) += M_ij;
			}
	}

	return density * m;
}

std::vector<int> Wedge::edge_to_node(int edge) {
	switch (edge) {
	case 0: return { 0, 1, 2 };
	case 1: return { 1, 2, 5, 4 };
	case 2: return { 2, 0, 3, 5 };
	case 3: return { 0, 1, 3, 4 };
	case 4: return { 3, 4, 5 };
	default:
		throw runtime_error("Error: wrong edge");
	}
}

std::array<double, 3> Wedge::normal(int edge) {
	std::vector<int> nodes = edge_to_node(edge);

	if (nodes.size() < 3) {
		throw std::runtime_error("Error: not enough nodes for normal calculation");
	}

	std::array<double, 3> p1 = { x[nodes[0]], y[nodes[0]], z[nodes[0]] };
	std::array<double, 3> p2 = { x[nodes[1]], y[nodes[1]], z[nodes[1]] };
	std::array<double, 3> p3 = { x[nodes[2]], y[nodes[2]], z[nodes[2]] };

	std::array<double, 3> v1 = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
	std::array<double, 3> v2 = { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] };

	std::array<double, 3> normal = {
		v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0]
	};

	return normal;
}

double Wedge::area_edge(int edge) {
	std::vector<int> nodes = edge_to_node(edge);

	if (nodes.size() < 3) {
		throw std::runtime_error("Error: not enough nodes for area calculation");
	}

	std::array<double, 3> p1 = { x[nodes[0]], y[nodes[0]], z[nodes[0]] };
	std::array<double, 3> p2 = { x[nodes[1]], y[nodes[1]], z[nodes[1]] };
	std::array<double, 3> p3 = { x[nodes[2]], y[nodes[2]], z[nodes[2]] };

	std::array<double, 3> v1 = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
	std::array<double, 3> v2 = { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] };

	std::array<double, 3> cross = {
		v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0]
	};

	double area = std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

	if (nodes.size() == 3)
		area /= 2.0;

	return area;
}

void Wedge::set_pressure(int edge, double value) {
	std::vector<int> nodes = edge_to_node(edge);

	std::array<double, 3> n = normal(edge);
	double len_norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

	for (auto& component : n)
		component /= len_norm;

	double area = area_edge(edge);

	for (int i = 0; i < 3; i++) {
		std::pair <int, int> pair(edge, i);
		load.insert({ pair, -value * n[i] * area });
	}
}

double Wedge::Volume() {
	std::complex<double> S = 0;

	for (int i = 0; i < 6; i++)
		for (int gp = 0; gp < 6; gp++) {
			double ksi = gaussPoint(KSI, gp);
			double eta = gaussPoint(ETA, gp);
			double zeta = gaussPoint(ZETA, gp);

			S += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
				FF(ksi, eta, zeta)[i] * std::abs(J(ksi, eta, zeta).determinant());
		}

	return S.real();
}

double Wedge::gaussPoint(LocVar var, int i) {
	std::vector<std::vector<double>> gp = {
		{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 },
		{ 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 },
		{ -0.577350269189626, -0.577350269189626, -0.577350269189626,
		   0.577350269189626,  0.577350269189626,  0.577350269189626 }
	};
	return gp[static_cast<int>(var)][i];
}

double Wedge::weight(LocVar var, int i) {
	double tri_w[3] = { 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 };
	// quad w = 1.0
	return tri_w[i % 3];
}

bool Wedge::pointInElem(std::vector<double> point) {
	return false;
}

std::vector<double> Wedge::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}