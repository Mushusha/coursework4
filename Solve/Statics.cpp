#include "Statics.h"

void Statics::calcDisp() {
	zeroDiagonalCheck();

    try {
        //Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;

        if (K.rows() != K.cols())
            throw runtime_error("Matrix is not square");

        if (K.nonZeros() == 0)
            throw runtime_error("Matrix is empty");

        solver.compute(K);

        if (solver.info() != Eigen::Success) {
            switch (solver.info()) {
            case Eigen::NumericalIssue:
                throw runtime_error("Numerical issue in decomposition");
            case Eigen::NoConvergence:
                throw runtime_error("No convergence in decomposition");
            case Eigen::InvalidInput:
                throw runtime_error("Invalid input matrix");
            default:
                throw runtime_error("Unknown error in decomposition");
            }
        }

        U = solver.solve(F);

        if (solver.info() != Eigen::Success) {
            throw runtime_error("Error in solve phase");
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Solver error: " << e.what() << std::endl;
    }
}

