#include "QuadElements.h"

std::vector<double> quadElement::FF(double ksi, double eta, double zeta) {
	std::vector<double> FF;
	FF.resize(4);
	FF[0] = (1 - ksi) * (1 - eta) / 4;
	FF[1] = (1 + ksi) * (1 - eta) / 4;
	FF[2] = (1 + ksi) * (1 + eta) / 4;
	FF[3] = (1 - ksi) * (1 + eta) / 4;
	return FF;
}

Eigen::MatrixXd quadElement::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXd gradFF = Eigen::MatrixXd::Zero(2, 4);
	double h = 0.01;
	for (int i = 0; i < 4; i++) {
		gradFF(KSI, i) = (FF(ksi + h, eta)[i] - FF(ksi - h, eta)[i]) / (2 * h);
		gradFF(ETA, i) = (FF(ksi, eta + h)[i] - FF(ksi, eta - h)[i]) / (2 * h);
	}
	return gradFF;
}

Eigen::MatrixXd quadElement::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, 2);
	for (int i = 0; i < 4; i++) {
		J(0, 0) += gradFF(ksi, eta)(KSI, i) * x[i];
		J(0, 1) += gradFF(ksi, eta)(KSI, i) * y[i];
		J(1, 0) += gradFF(ksi, eta)(ETA, i) * x[i];
		J(1, 1) += gradFF(ksi, eta)(ETA, i) * y[i];
	}
	return J;
}

Eigen::MatrixXd  quadElement::B(double ksi, double eta, double zeta) {
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

bool quadElement::pointInElem(std::vector<double> point) {
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

void quadElement::set_pressure(int edge, double value) {
	int comp = 0;
	if (edge == 0 || edge == 2)
		comp = 1;
	std::pair <int, int> pair(edge, comp);
	load.insert({ pair, -value });
}

double quadElement::len_edge(int edge) {
	std::pair<int, int> coords = edge_to_node(edge);
	return std::sqrt(std::pow((x[coords.first] - x[coords.second]), 2) + std::pow((y[coords.first] - y[coords.second]), 2));
}

std::pair<int, int> quadElement::edge_to_node(int edge) {
	if (edge != 3)
		return std::pair<int, int>(edge, edge + 1);
	else if (edge == 3)
		return std::pair<int, int>(edge, 0);
	else
		throw runtime_error("Error: wrong edge");
}

Eigen::MatrixXd quadElement::localK() {

	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(8, 8);
	std::vector <double> ksi = { -0.5774, 0.5774, -0.5774, -0.5774 };
	std::vector <double> eta = { -0.5774, -0.5774, -0.5774, 0.5774 };

	for (int gp = 0; gp < 4; gp++)
		for (int i = 0; i < 4; i++) {
			Eigen::MatrixXd B_i = Eigen::MatrixXd::Zero(3, 2);
			B_i(0, 0) = B(ksi[gp], eta[gp])(0, 2 * i);
			B_i(1, 1) = B(ksi[gp], eta[gp])(1, 2 * i + 1);
			B_i(2, 0) = B(ksi[gp], eta[gp])(2, 2 * i);
			B_i(2, 1) = B(ksi[gp], eta[gp])(2, 2 * i + 1);
			for (int j = 0; j < 4; j++) {
				Eigen::MatrixXd B_j = Eigen::MatrixXd::Zero(3, 2);
				B_j(0, 0) = B(ksi[gp], eta[gp])(0, 2 * i);
				B_j(1, 1) = B(ksi[gp], eta[gp])(1, 2 * i + 1);
				B_j(2, 0) = B(ksi[gp], eta[gp])(2, 2 * i);
				B_j(2, 1) = B(ksi[gp], eta[gp])(2, 2 * i + 1);

				Eigen::MatrixXd k_ij = B_i.transpose() * planeStrainD() * B_j * std::abs(J(eta[gp], ksi[gp]).determinant());

				k(2 * i, 2 * j) += k_ij(0, 0);
				k(2 * i, 2 * j + 1) += k_ij(0, 1);
				k(2 * i + 1, 2 * j) += k_ij(1, 0);
				k(2 * i + 1, 2 * j + 1) += k_ij(1, 1);

				//k.coeffRef(i, j) += B(ksi[gp], eta[gp]).transpose() * planeStrainD() * B(ksi[gp], eta[gp])* std::abs(J(eta[gp], ksi[gp]).determinant());
		}
	}
	return k;
}

std::vector<double> quadElement::localF() {
	std::vector<double> F;
	F.resize(8);

	// l.first.first - edge, l.first.second - comp, l.second - value
	for (auto const& l : load) {
		std::pair<int, int> node = edge_to_node(l.first.first);
		F[2 * node.first + l.first.second] += l.second * len_edge(l.first.first) / 2;
		F[2 * node.second + l.first.second] += l.second * len_edge(l.first.first) / 2;
	}
	return F;
}

Eigen::MatrixXd quadElement::localC() {
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(4, 4);
	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774 };

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int gp = 0; gp < 4; gp++)
//				std::cout << FF(ksi[gp], eta[gp])[i] * FF(ksi[gp], eta[gp])[j] * std::abs(J(eta[gp], ksi[gp]).determinant()) << std::endl;
				c(i, j) += FF(ksi[gp], eta[gp])[i] * FF(ksi[gp], eta[gp])[j] * std::abs(J(eta[gp], ksi[gp]).determinant());
	return c;
}

std::vector<double> quadElement::localR(std::vector<double> value) {
	std::vector<double> R;
	R.resize(4);

	std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774 };
	std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774 };
	for (int i = 0; i < 4; i++)
		for (int gp = 0; gp < 4; gp++)
			R[i] += value[i] * FF(ksi[gp], eta[gp])[i] * std::abs(J(eta[gp], ksi[gp]).determinant()); // ??

	return R;
}

//double quadInfN(double ksi, double eta, std::vector <double> a) {
//	return ((1 + eta) * (-ksi / (1 - ksi) * a[0] + (1 + ksi / (1 - ksi)) * a[1]) + (1 - eta) * (-ksi / (1 - ksi) * a[2] + (1 + ksi / (1 - ksi)) * a[3])) / 2;
//}
// 0 - x_C, 1 - x_Q, 2 - x_C1, 3 - x_Q1
