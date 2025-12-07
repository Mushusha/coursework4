#include "Hex.h"


std::vector<std::complex<double>> Hex::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
	FF.resize(8);
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
	Eigen::MatrixXcd gradFF = Eigen::MatrixXcd::Zero(3, 8);

	gradFF(KSI, 0) = -(1 - eta) * (1 - zeta) / 8;
	gradFF(KSI, 1) =  (1 - eta) * (1 - zeta) / 8;
	gradFF(KSI, 2) =  (1 + eta) * (1 - zeta) / 8;
	gradFF(KSI, 3) = -(1 + eta) * (1 - zeta) / 8;
	gradFF(KSI, 4) = -(1 - eta) * (1 + zeta) / 8;
	gradFF(KSI, 5) =  (1 - eta) * (1 + zeta) / 8;
	gradFF(KSI, 6) =  (1 + eta) * (1 + zeta) / 8;
	gradFF(KSI, 7) = -(1 + eta) * (1 + zeta) / 8;

	gradFF(ETA, 0) = -(1 - ksi) * (1 - zeta) / 8;
	gradFF(ETA, 1) = -(1 + ksi) * (1 - zeta) / 8;
	gradFF(ETA, 2) =  (1 + ksi) * (1 - zeta) / 8;
	gradFF(ETA, 3) =  (1 - ksi) * (1 - zeta) / 8;
	gradFF(ETA, 4) = -(1 - ksi) * (1 + zeta) / 8;
	gradFF(ETA, 5) = -(1 + ksi) * (1 + zeta) / 8;
	gradFF(ETA, 6) =  (1 + ksi) * (1 + zeta) / 8;
	gradFF(ETA, 7) =  (1 - ksi) * (1 + zeta) / 8;

	gradFF(ZETA, 0) = -(1 - ksi) * (1 - eta) / 8;
	gradFF(ZETA, 1) = -(1 + ksi) * (1 - eta) / 8;
	gradFF(ZETA, 2) = -(1 + ksi) * (1 + eta) / 8;
	gradFF(ZETA, 3) = -(1 - ksi) * (1 + eta) / 8;
	gradFF(ZETA, 4) =  (1 - ksi) * (1 - eta) / 8;
	gradFF(ZETA, 5) =  (1 + ksi) * (1 - eta) / 8;
	gradFF(ZETA, 6) =  (1 + ksi) * (1 + eta) / 8;
	gradFF(ZETA, 7) =  (1 - ksi) * (1 + eta) / 8;

	return gradFF;
}

Eigen::MatrixXcd Hex::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd J = Eigen::MatrixXcd::Zero(3, 3);
	Eigen::MatrixXcd grad = gradFF(ksi, eta, zeta);

	for (int i = 0; i < 8; i++) {
		J(0, 0) += grad(KSI, i) * x[i];
		J(0, 1) += grad(ETA, i) * x[i];
		J(0, 2) += grad(ZETA, i) * x[i];
		J(1, 0) += grad(KSI, i) * y[i];
		J(1, 1) += grad(ETA, i) * y[i];
		J(1, 2) += grad(ZETA, i) * y[i];
		J(2, 0) += grad(KSI, i) * z[i];
		J(2, 1) += grad(ETA, i) * z[i];
		J(2, 2) += grad(ZETA, i) * z[i];
	}
	return J;
}

Eigen::MatrixXcd Hex::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(6, 24);
	Eigen::Matrix3cd invJ;
	invJ = J(ksi, eta, zeta).inverse();
	Eigen::MatrixXcd dN = invJ.transpose() * gradFF(ksi, eta, zeta);

	for (int i = 0; i < 8; i++) {
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

Eigen::MatrixXcd Hex::localK() {
	Eigen::MatrixXcd k = Eigen::MatrixXcd::Zero(24, 24);

	for (int gp = 0; gp < 8; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
		double zeta = gaussPoint(ZETA, gp);

		k += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) * B(ksi, eta, zeta).transpose() * 
			D * B(ksi, eta, zeta) * std::abs(J(ksi, eta, zeta).determinant());
	}
	return k;
}

std::vector<double> Hex::localF(double mult) {
	std::vector<double> F;
	F.resize(24);

	// l.first.first - edge, l.first.second - comp, l.second - value
	if (load.size() != 0)
		for (auto const& l : load) {
			std::vector<int> node = edge_to_node(l.first.first);
			for (int i = 0; i < 4; i++)
				F[3 * node[i] + l.first.second] += mult * l.second / 4;
		}
	return F;
}

Eigen::MatrixXd Hex::localC() {
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(8, 8);

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			for (int gp = 0; gp < 8; gp++) {
				double ksi = gaussPoint(KSI, gp);
				double eta = gaussPoint(ETA, gp);
				double zeta = gaussPoint(ZETA, gp);

				c(i, j) += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) * FF(ksi, eta, zeta)[i].real() * 
					FF(ksi, eta, zeta)[j].real() * std::abs(J(ksi, eta, zeta).determinant());
			}
	return c;
}

std::vector<double> Hex::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(8);

	for (int i = 0; i < 8; i++)
		for (int gp = 0; gp < 8; gp++) {
			double ksi = gaussPoint(KSI, gp);
			double eta = gaussPoint(ETA, gp);
			double zeta = gaussPoint(ZETA, gp);

			R[i] += value[i] * weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) * 
				FF(ksi, eta, zeta)[i].real() * std::abs(J(ksi, eta, zeta).determinant());
		}
	return R;
}

Eigen::MatrixXcd Hex::localM() {
	if (density == 0.0)
		throw runtime_error("Error: density is zero in element " + to_string(id));

	Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(24, 24);

	for (size_t gp = 0; gp < 8; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
		double zeta = gaussPoint(ZETA, gp);

		for (int i = 0; i < 8; i++)
			for (int j = 0; j < 8; j++) {
				complex<double> M_ij = weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) * 
					FF(ksi, eta, zeta)[i] * FF(ksi, eta, zeta)[j] * std::abs(J(ksi, eta, zeta).determinant());

				m(3 * i, 3 * j) += M_ij;
				m(3 * i + 1, 3 * j + 1) += M_ij;
				m(3 * i + 2, 3 * j + 2) += M_ij;
			}
	}

	return density * m;
}

std::vector<int> Hex::edge_to_node(int edge) {
	switch (edge) {
	case 0:
		return { 1, 0, 3, 2 };
	case 1:
		return { 0, 1, 5, 4 };
	case 2:
		return { 1, 2, 6, 5 };
	case 3:
		return { 2, 3, 7, 6 };
	case 4:
		return { 3, 0, 4, 7 };
	case 5:
		return { 4, 5, 6, 7 };
	default:
		throw runtime_error("Error: wrong edge");
	}
}

std::array<double, 3> Hex::normal(int edge) {
	std::vector<int> nodes = edge_to_node(edge);

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

double Hex::area_edge(int edge) {
	std::vector<int> nodes = edge_to_node(edge);

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

	return area;
}

void Hex::set_pressure(int edge, double value) {
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

double Hex::Volume() {
	std::complex<double> S = 0;

	for (int i = 0; i < 8; i++)
		for (int gp = 0; gp < 8; gp++) {
			double ksi = gaussPoint(KSI, gp);
			double eta = gaussPoint(ETA, gp);
			double zeta = gaussPoint(ZETA, gp);

			S += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) * 
				FF(ksi, eta, zeta)[i] * std::abs(J(ksi, eta, zeta).determinant());
		}

	return S.real();
}

double Hex::gaussPoint(LocVar var, int i) {
	static const double gauss_coord = 1.0 / std::sqrt(3.0);
	
	static const std::array<std::array<double, 8>, 3> gp = {{
		{ -gauss_coord,  gauss_coord,  gauss_coord, -gauss_coord, -gauss_coord,  gauss_coord,  gauss_coord, -gauss_coord }, // KSI
		{ -gauss_coord, -gauss_coord,  gauss_coord,  gauss_coord, -gauss_coord, -gauss_coord,  gauss_coord,  gauss_coord }, // ETA
		{ -gauss_coord, -gauss_coord, -gauss_coord, -gauss_coord,  gauss_coord,  gauss_coord,  gauss_coord,  gauss_coord }  // ZETA
	}};

	return gp[static_cast<int>(var)][i];
}

double Hex::weight(LocVar var, int i) {
	return 1.0;
}

bool Hex::pointInElem(std::vector<double> point) {
	return false;
}

std::vector<double> Hex::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}