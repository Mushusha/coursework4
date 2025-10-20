#include "MathMV.h"

void print(Eigen::SparseMatrix<std::complex<double>>& A) {
	for (int i = 0; i < A.innerSize(); i++)
		for (int j = 0; j < A.outerSize(); j++)
			if (A.coeffRef(i, j) != 0.0)
				std::cout << i << " " << j << "\t" << A.coeffRef(i, j) << " " << std::endl;
}

void print(Eigen::SparseVector<std::complex<double>>& B) {
	for (int i = 0; i < B.size(); i++)
		if (B.coeffRef(i) != 0.0)
			std::cout << i << "\t" << B.coeffRef(i) << std::endl;
}

void print(Eigen::VectorX <std::complex<double>>& B) {
	for (int i = 0; i < B.size(); i++)
		if (B(i) != 0.0)
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
