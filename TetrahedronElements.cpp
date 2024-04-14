#include "TetrahedronElements.h"

Eigen::Matrix4f tetraMatrixC(std::vector <float> x, std::vector <float> y, std::vector <float> z) {

	Eigen::Matrix4f C;
	for (int i = 0; i < 4; i++) {
		C(i, 0) = 1;
		C(i, 1) = x[i];
		C(i, 2) = y[i];
		C(i, 3) = z[i];

	}
	return C;
}

Eigen::MatrixXf tetraMatrixB(std::vector <float> x, std::vector <float> y, std::vector <float> z) {

	Eigen::MatrixXf B = Eigen::MatrixXf::Zero(6, 12);
	Eigen::Matrix4f coef;
	Eigen::Matrix3f a, b, c, d;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			a(j, 0) = x[(i + j + 1) % 4];
			a(j, 1) = y[(i + j + 1) % 4];
			a(j, 2) = z[(i + j + 1) % 4];

			b(j, 0) = 1;
			b(j, 1) = y[(i + j + 1) % 4];
			b(j, 2) = z[(i + j + 1) % 4];

			c(j, 0) = x[(i + j + 1) % 4];
			c(j, 1) = 1;
			c(j, 2) = z[(i + j + 1) % 4];

			d(j, 0) = x[(i + j + 1) % 4];
			d(j, 1) = y[(i + j + 1) % 4];
			d(j, 2) = 1;
		}
		coef(0, i) = a.determinant();
		coef(1, i) = -1 * b.determinant();
		coef(2, i) = -1 * c.determinant();
		coef(3, i) = -1 * d.determinant();
	}

	for (int i = 0; i < 4; i++) {
		B(0, i * 3) = coef(i, 1);
		B(1, i * 3 + 1) = coef(i, 2);
		B(2, i * 3 + 2) = coef(i, 3);
		B(3, i * 3) = coef(i, 2);
		B(3, i * 3 + 1) = coef(i, 1);
		B(4, i * 3 + 1) = coef(i, 3);
		B(4, i * 3 + 2) = coef(i, 2);
		B(5, i * 3) = coef(i, 3);
		B(5, i * 3 + 2) = coef(i, 1);
	}

	return B / (6 * tetraMatrixC(x, y, z).determinant());
}

Eigen::MatrixXf tetraLocK(std::vector <float> x, std::vector <float> y, std::vector <float> z) {

	return tetraMatrixB(x, y, z).transpose() * threeMatrixD() * tetraMatrixB(x, y, z) * tetraMatrixC(x, y, z).determinant() / 6;

}