#include "MathMV.h"

void writeMatrix(Eigen::SparseMatrix<double> A) {
	for (int i = 0; i < A.innerSize(); i++) {
		for (int j = 0; j < A.outerSize(); j++)
			std::cout << A.coeffRef(i, j) << " ";
		std::cout << std::endl;
	}
}

void writeVector(Eigen::SparseVector<double> B) {
	for (int i = 0; i < B.size(); i++)
		std::cout << B.coeffRef(i) << std::endl;
}

void productMV(Eigen::MatrixXd A, std::vector<double> b, std::vector<double>& c) {
	if (A.outerSize() != b.size()) {
		std::cout << "Error: incorrect size" << std::endl;
		return;
	}
	c.resize(A.innerSize(), 0); 
	for (int row = 0; row < A.innerSize(); row++)
		for (int col = 0; col < A.outerSize(); col++)
			c[row] += b[col] * A(row, col);
}
