#include "Pyr.h"


std::vector<std::complex<double>> Pyr::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
	FF.resize(5);
	
	FF[0] = (1 - ksi) * (1 - eta) * (1 - zeta) / 8;
	FF[1] = (1 + ksi) * (1 - eta) * (1 - zeta) / 8;
	FF[2] = (1 + ksi) * (1 + eta) * (1 - zeta) / 8;
	FF[3] = (1 - ksi) * (1 + eta) * (1 - zeta) / 8;
	FF[4] = (1 + zeta) / 2;
	
	return FF;
}

Eigen::MatrixXcd Pyr::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd gradFF = Eigen::MatrixXcd::Zero(3, 5);
	
	gradFF(KSI, 0) = -(1 - eta) * (1 - zeta) / 8;
	gradFF(KSI, 1) =  (1 - eta) * (1 - zeta) / 8;
	gradFF(KSI, 2) =  (1 + eta) * (1 - zeta) / 8;
	gradFF(KSI, 3) = -(1 + eta) * (1 - zeta) / 8;
	gradFF(KSI, 4) =  0;
	
	gradFF(ETA, 0) = -(1 - ksi) * (1 - zeta) / 8;
	gradFF(ETA, 1) = -(1 + ksi) * (1 - zeta) / 8;
	gradFF(ETA, 2) =  (1 + ksi) * (1 - zeta) / 8;
	gradFF(ETA, 3) =  (1 - ksi) * (1 - zeta) / 8;
	gradFF(ETA, 4) =  0;
	
	gradFF(ZETA, 0) = -(1 - ksi) * (1 - eta) / 8;
	gradFF(ZETA, 1) = -(1 + ksi) * (1 - eta) / 8;
	gradFF(ZETA, 2) = -(1 + ksi) * (1 + eta) / 8;
	gradFF(ZETA, 3) = -(1 - ksi) * (1 + eta) / 8;
	gradFF(ZETA, 4) =  0.5;
	
	return gradFF;
}

Eigen::MatrixXcd Pyr::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd J = Eigen::MatrixXcd::Zero(3, 3);
	Eigen::MatrixXcd grad = gradFF(ksi, eta, zeta);
	
	for (int i = 0; i < 5; i++) {
		J(0, 0) += grad(KSI, i) * x[i];
		J(0, 1) += grad(KSI, i) * y[i];
		J(0, 2) += grad(KSI, i) * z[i];
		J(1, 0) += grad(ETA, i) * x[i];
		J(1, 1) += grad(ETA, i) * y[i];
		J(1, 2) += grad(ETA, i) * z[i];
		J(2, 0) += grad(ZETA, i) * x[i];
		J(2, 1) += grad(ZETA, i) * y[i];
		J(2, 2) += grad(ZETA, i) * z[i];
	}
	return J;
}

Eigen::MatrixXcd Pyr::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(6, 15);
	Eigen::Matrix3cd invJ;
	invJ = J(ksi, eta, zeta).inverse();
	Eigen::MatrixXcd dN = invJ * gradFF(ksi, eta, zeta);
	
	for (int i = 0; i < 5; i++) {
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

Eigen::MatrixXcd Pyr::localK() {
	Eigen::MatrixXcd k = Eigen::MatrixXcd::Zero(15, 15);
	
	for (int gp = 0; gp < 5; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
		double zeta = gaussPoint(ZETA, gp);
		
		k += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
			B(ksi, eta, zeta).transpose() * D * B(ksi, eta, zeta) *
			std::abs(J(ksi, eta, zeta).determinant());
	}
	return k;
}

std::vector<double> Pyr::localF(double mult) {
	std::vector<double> F;
	F.resize(15);
	
	// l.first.first - edge, l.first.second - comp, l.second - value
	if (load.size() != 0)
		for (auto const& l : load) {
			std::vector<int> node = edge_to_node(l.first.first);
			for (size_t i = 0; i < node.size(); i++)
				F[3 * node[i] + l.first.second] += mult * l.second / node.size();
		}
	return F;
}

Eigen::MatrixXd Pyr::localC() {
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(5, 5);
	
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			for (int gp = 0; gp < 5; gp++) {
				double ksi = gaussPoint(KSI, gp);
				double eta = gaussPoint(ETA, gp);
				double zeta = gaussPoint(ZETA, gp);
				
				c(i, j) += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
					FF(ksi, eta, zeta)[i].real() * FF(ksi, eta, zeta)[j].real() *
					std::abs(J(ksi, eta, zeta).determinant());
			}
	return c;
}

std::vector<double> Pyr::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(5);
	
	for (int i = 0; i < 5; i++)
		for (int gp = 0; gp < 5; gp++) {
			double ksi = gaussPoint(KSI, gp);
			double eta = gaussPoint(ETA, gp);
			double zeta = gaussPoint(ZETA, gp);
			
			R[i] += value[i] * weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
				FF(ksi, eta, zeta)[i].real() * std::abs(J(ksi, eta, zeta).determinant());
		}
	return R;
}

Eigen::MatrixXcd Pyr::localM() {
	if (density == 0.0)
		throw runtime_error("Error: density is zero in element " + to_string(id));
	
	Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(15, 15);
	
	for (size_t gp = 0; gp < 5; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
		double zeta = gaussPoint(ZETA, gp);
		
		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++) {
				complex<double> M_ij = weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
					FF(ksi, eta, zeta)[i] * FF(ksi, eta, zeta)[j] * std::abs(J(ksi, eta, zeta).determinant());
				
				m(3 * i, 3 * j) += M_ij;
				m(3 * i + 1, 3 * j + 1) += M_ij;
				m(3 * i + 2, 3 * j + 2) += M_ij;
			}
	}
	
	return density * m;
}

std::vector<int> Pyr::edge_to_node(int edge) {
	switch (edge) {
	case 0: return { 0, 1, 2, 3 };
	case 1: return { 0, 1, 4 };
	case 2: return { 1, 2, 4 };
	case 3: return { 2, 3, 4 };
	case 4: return { 3, 0, 4 };
	default:
		throw runtime_error("Error: wrong edge");
	}
}

std::array<double, 3> Pyr::normal(int edge) {
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

double Pyr::area_edge(int edge) {
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

void Pyr::set_pressure(int edge, double value) {
	std::vector<int> nodes = edge_to_node(edge);
	
	std::array<double, 3> n = normal(edge);
	double len_norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	
	for (auto& component : n)
		component /= len_norm;
	
	double area = area_edge(edge);
	
	for (int i = 0; i < 3; i++) {
		std::pair<int, int> pair(edge, i);
		load.insert({ pair, -value * n[i] * area });
	}
}

double Pyr::Volume() {
	std::complex<double> S = 0;
	
	for (int gp = 0; gp < 5; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
		double zeta = gaussPoint(ZETA, gp);
		
		S += weight(KSI, gp) * weight(ETA, gp) * weight(ZETA, gp) *
			std::abs(J(ksi, eta, zeta).determinant());
	}
	
	return S.real();
}

double Pyr::gaussPoint(LocVar var, int i) {

	static const double a = 0.584237394672177;
	static const double z0 = -2.0 / 3.0;
	static const double z1 = 1.0 / 5.0;
	
	std::vector<std::vector<double>> gp = {
		{ 0.0,  -a,   a,   a,  -a },
		{ 0.0,  -a,  -a,   a,   a },
		{ z0,   z1,  z1,  z1,  z1 }
	};
	
	return gp[static_cast<int>(var)][i];
}

double Pyr::weight(LocVar var, int i) {
	static const double w0 = 81.0 / 200.0; 
	static const double w1 = 125.0 / 27.0 / 8.0;
	
	static const double weights[] = { w0, w1, w1, w1, w1 };
	return weights[i];
}

bool Pyr::pointInElem(std::vector<double> point) {
	return false;
}

std::vector<double> Pyr::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}

