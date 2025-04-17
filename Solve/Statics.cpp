#include "Statics.h"

void Statics::calcDisp() {
	zeroDiagonalCheck();

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
	solver.compute(K);

	if (solver.info() != Eigen::Success)
		throw runtime_error("Error in K");

	U = solver.solve(F);
}

