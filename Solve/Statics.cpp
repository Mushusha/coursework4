#include "Statics.h"

void Statics::calcDisp() {
	zeroDiagonalCheck();
	printMatrix(K);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(K);

	if (solver.info() != Eigen::Success)
		throw runtime_error("Error in K");

	U = solver.solve(F);
}

