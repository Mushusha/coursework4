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
	double h = 0.01;
	for (int i = 0; i < 8; i++) {
		gradFF(KSI, i) = (FF(ksi + h, eta, zeta)[i] - FF(ksi - h, eta, zeta)[i]) / (2 * h);
		gradFF(ETA, i) = (FF(ksi, eta + h, zeta)[i] - FF(ksi, eta - h, zeta)[i]) / (2 * h);
		gradFF(ZETA, i) = (FF(ksi, eta, zeta + h)[i] - FF(ksi, eta, zeta - h)[i]) / (2 * h);
	}
	return gradFF;
}

Eigen::MatrixXcd Hex::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd J = Eigen::MatrixXcd::Zero(3, 3);

	for (int i = 0; i < 8; i++) {
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
			F[2 * node[0] + l.first.second] += mult * l.second / 4;
			F[2 * node[1] + l.first.second] += mult * l.second / 4;
			F[2 * node[2] + l.first.second] += mult * l.second / 4;
			F[2 * node[3] + l.first.second] += mult * l.second / 4;
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

void Hex::set_pressure(int edge, double value) {
	//std::pair<int, int> node = edge_to_node(edge);
	//std::array<double, 2> comp;
	//comp[0] = -y[node.first] + y[node.second];
	//comp[1] = x[node.first] - x[node.second];
	//
	//if ((x[node.first] - x[node.second]) * (y[node.first] - y[(edge + 2) % 4]) -
	//	(y[node.first] - y[node.second]) * (x[node.first] - x[(edge + 2) % 4]) < 0)
	//	for (auto& i : comp)
	//		i *= -1;
	//
	//for (int i = 0; i < 2; i++) {
	//	std::pair <int, int> pair(edge, i);
	//	load.insert({ pair, -value * comp[i] / len_edge(edge) });
	//}
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
	std::vector<std::vector<double>> gp =
	{ { -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926 },
	{ -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926 },
	{ -0.57735026918926, -0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, 0.57735026918926, 0.57735026918926 } };

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