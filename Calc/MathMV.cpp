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
