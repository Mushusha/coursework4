#include "MathMV.h"
#include <limits>
#include <cmath>

void print(Eigen::SparseMatrix<std::complex<double>>& A) {
	for (int i = 0; i < A.innerSize(); i++)
		for (int j = 0; j < A.outerSize(); j++)
			if ((A.coeffRef(i, j).real() > 1e-10) || (A.coeffRef(i, j).real() < -1e-10))
				std::cout << i << " " << j << "\t" << A.coeffRef(i, j) << " " << std::endl;
}

void print(Eigen::SparseVector<std::complex<double>>& B) {
	for (int i = 0; i < B.size(); i++)
		if ((B.coeffRef(i).real() > 1e-10) || (B.coeffRef(i).real() < -1e-10))
			std::cout << i << "\t" << B.coeffRef(i) << std::endl;
}

void print(Eigen::VectorX <std::complex<double>>& B) {
	for (int i = 0; i < B.size(); i++)
		if ((B(i).real() > 1e-10) || (B(i).real() < -1e-10))
			std::cout << i << "\t" << B(i) << std::endl;
}

double line(std::vector<double> A, std::vector<double> B, std::vector<double> X) {
	// return value y - f(x)
	return (X[1] - A[1]) * (B[0] - A[0]) - (X[0] - A[0]) * (B[1] - A[1]);
}

double berlage(double omega, double A, double time) {
	double pi = 3.14159265358979323;
	double omega0 = 2 * pi * omega;
	double omega1 = omega0 / sqrt(3);
	double mult = A * pow(omega1, 2) * exp(-omega1 * time) / 4;
	double term1 = sin(omega0 * time) * (-1 * pow(time, 2) / omega1 + time / pow(omega1, 2) + 1 / pow(omega1, 3));
	double term2 = cos(omega0 * time) * (pow(time, 2) / omega1 + time / pow(omega1, 2));

	return mult * (term1 - sqrt(3) * term2);
}

double min_edge_length(const std::shared_ptr<Element>& elem) {
	int nodes_count = elem->nodes_count();
	double min_len = std::numeric_limits<double>::max();
	
	for (int i = 0; i < nodes_count; i++) {
		for (int j = i + 1; j < nodes_count; j++) {
			double dx = elem->get_coord(i, 0) - elem->get_coord(j, 0);
			double dy = elem->get_coord(i, 1) - elem->get_coord(j, 1);
			double dz = elem->get_coord(i, 2) - elem->get_coord(j, 2);
			
			double len = std::sqrt(dx * dx + dy * dy + dz * dz);
			if (len < min_len) {
				min_len = len;
			}
		}
	}
	
	return min_len;
}

void compute_gll_nodes_weights(int order, std::vector<double>& points, std::vector<double>& weights) {
	if (order < 1)
		throw std::runtime_error("Order must be at least 1");

	int n = order + 1;
	points.resize(n);
	weights.resize(n);

	points[0] = -1.0;
	points[n - 1] = 1.0;

	if (order == 1) {
		weights[0] = 1.0;
		weights[1] = 1.0;
		return;
	}

	if (order > 1) {
		for (int i = 1; i < n - 1; i++) {
			double x = -cos(3.14159265358979323846 * i / order);
			double delta;
			int iter = 0, max_iter = 100;

			do {
				double P_prev = 1.0;
				double P_curr = x;
				double P_next;

				for (int k = 2; k <= order; k++) {
					P_next = ((2.0 * k - 1.0) * x * P_curr - (k - 1.0) * P_prev) / k;
					P_prev = P_curr;
					P_curr = P_next;
				}

				double derivative = order * (P_prev - x * P_curr) / (1.0 - x * x);
				double second_deriv = (2.0 * x * derivative - order * (order + 1) * P_curr) / (1.0 - x * x);

				delta = -derivative / second_deriv;
				x += delta;
				iter++;
			} while (std::abs(delta) > 1e-14 && iter < max_iter);

			points[i] = x;
		}
	}

	std::sort(points.begin(), points.end());

	for (int i = 0; i < n; i++) {
		double x = points[i];

		double P_prev = 1.0;
		double P_curr = x;
		double P_next;

		for (int k = 2; k <= order; k++) {
			P_next = ((2.0 * k - 1.0) * x * P_curr - (k - 1.0) * P_prev) / k;
			P_prev = P_curr;
			P_curr = P_next;
		}

		weights[i] = 2.0 / (order * (order + 1) * P_curr * P_curr);
	}
}

void compute_gauss_radau_nodes_weights(int n, std::vector<double>& points, std::vector<double>& weights) {
	
	if (n < 1)
		throw std::runtime_error("Gauss-Radau order must be at least 1");

	points.resize(n);
	weights.resize(n);

	points[0] = -1.0;

	if (n == 1) {
		weights[0] = 2.0;
		return;
	}

	if (n == 2) {
		points[0] = -1.0;
		points[1] = 1.0 / 3.0;
		weights[0] = 0.5;
		weights[1] = 1.5;
		return;
	}

	if (n == 3) {
		points[0] = -1.0;
		points[1] = (1.0 - std::sqrt(6.0)) / 5.0;
		points[2] = (1.0 + std::sqrt(6.0)) / 5.0;
		weights[0] = 2.0 / 9.0;
		weights[1] = (16.0 + std::sqrt(6.0)) / 18.0;
		weights[2] = (16.0 - std::sqrt(6.0)) / 18.0;
		return;
	}

	if (n == 4) {
		points[0] = -1.0;
		points[1] = -0.575318923521694;
		points[2] = 0.181066271118531;
		points[3] = 0.822824080974592;
		weights[0] = 0.125;
		weights[1] = 0.657688639960120;
		weights[2] = 0.776386937686343;
		weights[3] = 0.440924422353537;
		return;
	}

	if (n == 5) {
		points[0] = -1.0;
		points[1] = -0.720480271312439;
		points[2] = -0.167180864737834;
		points[3] = 0.446313972723752;
		points[4] = 0.885791607770964;
		weights[0] = 0.08;
		weights[1] = 0.446207802167141;
		weights[2] = 0.623653045951483;
		weights[3] = 0.562712030298924;
		weights[4] = 0.287427121682452;
		return;
	}

	for (int i = 1; i < n; i++) {
		double x = -std::cos(3.14159265358979323846 * (2.0 * i - 1.0) / (2.0 * (n - 1)));
		double delta;
		int iter = 0, max_iter = 100;

		do {
			double P_prev = 1.0;
			double P_curr = x;
			double P_next;

			for (int k = 2; k <= n; k++) {
				P_next = ((2.0 * k - 1.0) * x * P_curr - (k - 1.0) * P_prev) / k;
				P_prev = P_curr;
				P_curr = P_next;
			}

			double Pn = P_curr;
			double Pn_1 = P_prev;

			double f = (1.0 + x) * Pn_1 + Pn;

			double dPn_1 = (n - 1) * (x * Pn_1 - P_prev) / (x * x - 1.0);
			if (std::abs(x * x - 1.0) < 1e-10) {
				dPn_1 = (n - 1) * n * Pn_1 / 2.0;
			}
			double dPn = n * (x * Pn - Pn_1) / (x * x - 1.0);
			if (std::abs(x * x - 1.0) < 1e-10) {
				dPn = n * (n + 1) * Pn / 2.0;
			}

			double df = Pn_1 + (1.0 + x) * dPn_1 + dPn;

			delta = -f / df;
			x += delta;
			iter++;
		} while (std::abs(delta) > 1e-14 && iter < max_iter);

		points[i] = x;
	}

	std::sort(points.begin(), points.end());

	for (int i = 0; i < n; i++) {
		double x = points[i];

		double P_prev = 1.0;
		double P_curr = x;
		double P_next;

		for (int k = 2; k <= n; k++) {
			P_next = ((2.0 * k - 1.0) * x * P_curr - (k - 1.0) * P_prev) / k;
			P_prev = P_curr;
			P_curr = P_next;
		}

		double Pn = P_curr;

		if (i == 0) {
			weights[i] = 2.0 / (n * n);
		}
		else {
			weights[i] = (1.0 - x) / (n * n * Pn * Pn);
		}
	}
}