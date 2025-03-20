#include "Quad.h"

std::vector<double> Quad::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF.resize(4);
	FF[0] = (1 - ksi) * (1 - eta) / 4;
	FF[1] = (1 + ksi) * (1 - eta) / 4;
	FF[2] = (1 + ksi) * (1 + eta) / 4;
	FF[3] = (1 - ksi) * (1 + eta) / 4;
	return FF;
}

Eigen::MatrixXd Quad::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXd gradFF = Eigen::MatrixXd::Zero(2, 4);
	double h = 0.01;
	for (int i = 0; i < 4; i++) {
		gradFF(KSI, i) = (FF(ksi + h, eta)[i] - FF(ksi - h, eta)[i]) / (2 * h);
		gradFF(ETA, i) = (FF(ksi, eta + h)[i] - FF(ksi, eta - h)[i]) / (2 * h);
	}
	return gradFF;
}

Eigen::MatrixXd Quad::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, 2);

	for (int i = 0; i < 4; i++) {
		J(0, 0) += gradFF(ksi, eta)(KSI, i) * x[i];
		J(0, 1) += gradFF(ksi, eta)(KSI, i) * y[i];
		J(1, 0) += gradFF(ksi, eta)(ETA, i) * x[i];
		J(1, 1) += gradFF(ksi, eta)(ETA, i) * y[i];
	}
	return J;
}

Eigen::MatrixXd  Quad::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
	Eigen::Matrix2d invJ;
	invJ = J(ksi, eta).inverse();
	Eigen::MatrixXd dN = invJ * gradFF(ksi, eta);

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

Eigen::MatrixXd Quad::localK() {

	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(8, 8);
	std::vector <double> ksi = { -0.5774, 0.5774, 0.5774, -0.5774 };
	std::vector <double> eta = { -0.5774, -0.5774, 0.5774, 0.5774 };

	for (int gp = 0; gp < 4; gp++)
		k += B(ksi[gp], eta[gp]).transpose() * D * B(ksi[gp], eta[gp]) * std::abs(J(eta[gp], ksi[gp]).determinant());

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
	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int gp = 0; gp < 4; gp++)
				c(i, j) += FF(ksi[gp], eta[gp])[i] * FF(ksi[gp], eta[gp])[j] * std::abs(J(eta[gp], ksi[gp]).determinant());
	return c;
}

std::vector<double> Quad::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(4);

	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774 };
	for (int i = 0; i < 4; i++)
		for (int gp = 0; gp < 4; gp++)
			R[i] += value[i] * FF(ksi[gp], eta[gp])[i] * std::abs(J(eta[gp], ksi[gp]).determinant());
	return R;
}

Eigen::MatrixXd Quad::localM() {
	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(8, 8);
	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			if ((i + j) % 2 == 1)
				m(i, j) = 0;
			else
				for (int gp = 0; gp < 4; gp++) {
					m(2 * i, 2 * j) += FF(ksi[gp], eta[gp])[i] * FF(ksi[gp], eta[gp])[j] * std::abs(J(eta[gp], ksi[gp]).determinant());
					m(2 * i + 1, 2 * j + 1) += m(2 * i, 2 * j);
				}
	return density * m;
}

std::vector<double> Quad::coordFF(double x0, double y0, double z0) {
	Eigen::Vector2d F;
	std::vector x1 = { 1, 0, 3, 4 };
	std::vector y1 = { 0, 2, 4, 0 };

	double ksi0 = 0, eta0 = 0;
	for (int i = 0; i < 50; i++) {
		F(0) = 0;
		F(1) = 0;
		for (int j = 0; j < 4; j++) {
			F(0) += FF(ksi0, eta0)[j] * x[j];
			F(1) += FF(ksi0, eta0)[j] * y[j];
		}

		F(0) -= x0;
		F(1) -= y0;

		Eigen::Vector2d delta;
		delta = -1 * (J(ksi0, eta0).transpose()).inverse() * F;

		ksi0 += delta(0);
		eta0 += delta(1);

		if (std::hypot(delta(0), delta(1)) < 1e-8 || std::hypot(F(0), F(1)) < 1e-8)
			break;
	}

	std::vector<double> coord = { ksi0, eta0 };
	return coord;
}

double Quad::Volume() {
	double S = 0;
	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < 4; i++)
		for (int gp = 0; gp < 4; gp++)
			S += FF(ksi[gp], eta[gp])[i] * std::abs(J(eta[gp], ksi[gp]).determinant());

	return S;
}

//double quadInfN(double ksi, double eta, std::vector <double> a) {
//	return ((1 + eta) * (-ksi / (1 - ksi) * a[0] + (1 + ksi / (1 - ksi)) * a[1]) + (1 - eta) * (-ksi / (1 - ksi) * a[2] + (1 + ksi / (1 - ksi)) * a[3])) / 2;
//}
// 0 - x_C, 1 - x_Q, 2 - x_C1, 3 - x_Q1
