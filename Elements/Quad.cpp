#include "Quad.h"

std::vector<std::complex<double>> Quad::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> FF;
	FF.resize(4);
	FF[0] = (1 - ksi) * (1 - eta) / 4;
	FF[1] = (1 + ksi) * (1 - eta) / 4;
	FF[2] = (1 + ksi) * (1 + eta) / 4;
	FF[3] = (1 - ksi) * (1 + eta) / 4;
	return FF;
}

Eigen::MatrixXcd Quad::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd gradFF = Eigen::MatrixXcd::Zero(2, 4);
	double h = 0.01;
	for (int i = 0; i < 4; i++) {
		gradFF(KSI, i) = (FF(ksi + h, eta)[i] - FF(ksi - h, eta)[i]) / (2 * h);
		gradFF(ETA, i) = (FF(ksi, eta + h)[i] - FF(ksi, eta - h)[i]) / (2 * h);
	}
	return gradFF;
}

Eigen::MatrixXcd Quad::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd J = Eigen::MatrixXcd::Zero(2, 2);
	
	for (int i = 0; i < 4; i++) {
		J(0, 0) += gradFF(ksi, eta)(KSI, i) * x[i];
		J(0, 1) += gradFF(ksi, eta)(KSI, i) * y[i];
		J(1, 0) += gradFF(ksi, eta)(ETA, i) * x[i];
		J(1, 1) += gradFF(ksi, eta)(ETA, i) * y[i];
	}
	return J;
}

double Quad::gaussPoint(LocVar var, int i) {
	std::vector<std::vector<double>> gp = { { -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926 },
										  { -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926 },
										  { 0.0, 0.0, 0.0, 0.0 } };

	return gp[static_cast<int>(var)][i];
}

double Quad::weight(LocVar var, int i) {
	return 1.0;
}

Eigen::MatrixXcd Quad::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(3, 8);
	Eigen::Matrix2cd invJ;
	invJ = J(ksi, eta).inverse();
	Eigen::MatrixXcd dN = invJ * gradFF(ksi, eta);
	
	for (int i = 0; i < 4; i++) {
		B(0, 2 * i) = dN(X, i);
		B(1, 2 * i + 1) = dN(Y, i);
		B(2, 2 * i) = dN(Y, i);
		B(2, 2 * i + 1) = dN(X, i);
	}
	return B;
}

bool Quad::pointInElem(std::vector<double> point) {
	bool answer = true;
	for (int i = 0; i < 4; i++) {
		std::vector<double> a{ x[i], y[i] };
		std::vector<double> b{ x[(i + 1) % 4], y[(i + 1) % 4] };
		std::vector<double> c{ x[(i + 2) % 4], y[(i + 2) % 4] };
		std::vector<double> d{ x[(i + 3) % 4], y[(i + 3) % 4] };

		if (a == point || b == point || c == point || d == point)
			return true;
		if ((line(a, b, point) >= 0 && line(b, c, point) >= 0) &&
			(line(c, d, point) >= 0 && line(d, a, point) >= 0))
			answer &= true;
		else if ((line(a, b, point) <= 0 && line(b, c, point) <= 0) &&
			(line(c, d, point) <= 0 && line(d, a, point) <= 0))
			answer &= true;
		else
			return false;
	}
	return answer;
}

void Quad::set_pressure(int edge, double value) {
	std::pair<int, int> node = edge_to_node(edge);
	std::array<double, 2> comp;
	comp[0] = -y[node.first] + y[node.second];
	comp[1] = x[node.first] - x[node.second];

	if ((x[node.first] - x[node.second]) * (y[node.first] - y[(edge + 2) % 4]) -
		(y[node.first] - y[node.second]) * (x[node.first] - x[(edge + 2) % 4]) < 0)
		for (auto& i : comp)
			i *= -1;

	for (int i = 0; i < 2; i++) {
		std::pair <int, int> pair(edge, i);
		load.insert({ pair, -value * comp[i] / len_edge(edge) });
	}
}

double Quad::len_edge(int edge) {
	std::pair<int, int> coords = edge_to_node(edge);
	return std::sqrt(std::pow((x[coords.first] - x[coords.second]), 2) + std::pow((y[coords.first] - y[coords.second]), 2));
}

std::pair<int, int> Quad::edge_to_node(int edge) {
	if (edge != 3)
		return std::pair<int, int>(edge, edge + 1);
	else if (edge == 3)
		return std::pair<int, int>(edge, 0);
	else
		throw runtime_error("Error: wrong edge");
}

Eigen::MatrixXcd Quad::localK() {

	Eigen::MatrixXcd k = Eigen::MatrixXcd::Zero(8, 8);
	
	for (int gp = 0; gp < 4; gp++) {
		double ksi = gaussPoint(KSI, gp);
		double eta = gaussPoint(ETA, gp);
	
		k += weight(KSI, gp) * weight(ETA, gp) * B(ksi, eta).transpose() * D * B(ksi, eta) * std::abs(J(ksi, eta).determinant());
	}
	return k;
}

std::vector<double> Quad::localF() {
	std::vector<double> F;
	F.resize(8);

	// l.first.first - edge, l.first.second - comp, l.second - value
	if (load.size() != 0)
		for (auto const& l : load) {
			std::pair<int, int> node = edge_to_node(l.first.first);
			F[2 * node.first + l.first.second] += l.second * len_edge(l.first.first) / 2;
			F[2 * node.second + l.first.second] += l.second * len_edge(l.first.first) / 2;
		}
	return F;
}

Eigen::MatrixXd Quad::localC() {
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(4, 4);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int gp = 0; gp < 4; gp++) {
				double ksi = gaussPoint(KSI, gp);
				double eta = gaussPoint(ETA, gp);

				c(i, j) += weight(KSI, gp) * weight(ETA, gp) * FF(ksi, eta)[i].real() * FF(ksi, eta)[j].real() * std::abs(J(ksi, eta).determinant());
			}
	return c;
}

std::vector<double> Quad::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(4);

	for (int i = 0; i < 4; i++)
		for (int gp = 0; gp < 4; gp++) {
			double ksi = gaussPoint(KSI, gp);
			double eta = gaussPoint(ETA, gp);

			R[i] += value[i] * weight(KSI, gp) * weight(ETA, gp) * FF(ksi, eta)[i].real() * std::abs(J(ksi, eta).determinant());
		}
	return R;
}

Eigen::MatrixXcd Quad::localM() {
	if (density == 0.0)
		throw runtime_error("Error: density is zero in element " + to_string(id));

	Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(8, 8);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			if ((i + j) % 2 == 1)
				m(i, j) = 0;
			else
				for (int gp = 0; gp < 4; gp++) {
					double ksi = gaussPoint(KSI, gp);
					double eta = gaussPoint(ETA, gp);

					m(2 * i, 2 * j) += weight(KSI, gp) * weight(ETA, gp) * FF(ksi, eta)[i] * FF(ksi, eta)[j] * std::abs(J(ksi, eta).determinant());
					
					m(2 * i + 1, 2 * j + 1) += m(2 * i, 2 * j);
				}
	return density * m;
}

std::vector<double> Quad::coordFF(double x0, double y0, double z0) {
	Eigen::Vector2d F;
	double ksi0 = 0, eta0 = 0;
	for (int i = 0; i < 50; i++) {
		F(0) = 0;
		F(1) = 0;
		for (int j = 0; j < 4; j++) {
			F(0) += FF(ksi0, eta0)[j].real() * x[j];
			F(1) += FF(ksi0, eta0)[j].real() * y[j];
		}

		F(0) -= x0;
		F(1) -= y0;

		Eigen::Vector2cd delta;
		delta = -1 * (J(ksi0, eta0).transpose()).inverse() * F;

		ksi0 += delta(0).real();
		eta0 += delta(1).real();

		if (std::hypot(delta(0).real(), delta(1).real()) < 1e-8 || std::hypot(F(0), F(1)) < 1e-8)
			break;
	}

	std::vector<double> coord = { ksi0, eta0 };
	return coord;
}

double Quad::Volume() {
	std::complex<double> S = 0;
	
	for (int i = 0; i < 4; i++)
		for (int gp = 0; gp < 4; gp++) {
			double ksi = gaussPoint(KSI, gp);
			double eta = gaussPoint(ETA, gp);
	
			S += weight(KSI, gp) * weight(ETA, gp) * FF(ksi, eta)[i] * std::abs(J(ksi, eta).determinant());
		}
	
	return S.real();
}