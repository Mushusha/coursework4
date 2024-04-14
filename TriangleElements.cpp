#include "TriangleElements.h"

triangleElement::triangleElement(std::vector <int> _nodes, std::vector <float> _x, std::vector <float> _y) {
	for (int i = 0; i < _nodes.size(); i++) {
		nodes.push_back(_nodes[i]);
		x.push_back(_x[i]);
		y.push_back(_y[i]);
		z.push_back(0);
	}
	locK = triangleLocK(x, y);
}

Eigen::Matrix3f triangleMatrixC(std::vector <float> x, std::vector <float> y) {

	Eigen::Matrix3f C;
	for (int i = 0; i < 3; i++) {
		C(i, 0) = 1;
		C(i, 1) = x[i];
		C(i, 2) = y[i];
	}
	return C;
}

Eigen::MatrixXf triangleMatrixB(std::vector <float> x, std::vector <float> y) {

	Eigen::MatrixXf B(3, 6);
	for (int i = 0; i < 3; i++) {
		B(0, 2 * i) = y[(1 + i) % 3] - y[(2 + i) % 3];
		B(0, 2 * i + 1) = 0;
		B(1, 2 * i) = 0;
		B(1, 2 * i + 1) = x[(2 + i) % 3] - x[(1 + i) % 3];
		B(2, 2 * i) = x[(2 + i) % 3] - x[(1 + i) % 3];
		B(2, 2 * i + 1) = y[(1 + i) % 3] - y[(2 + i) % 3];
	}
	B = B / triangleMatrixC(x, y).determinant() / 2;
	return B;
}

Eigen::MatrixXf triangleLocK(std::vector <float> x, std::vector <float> y) {

	return triangleMatrixB(x, y).transpose() * twoMatrixD() * triangleMatrixB(x, y) * triangleMatrixC(x, y).determinant() / 2;
}