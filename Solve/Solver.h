#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/Core"
#include "unsupported/Eigen/IterativeSolvers"

#include "Calc/Data.h"

class Solver {
public:
	Solver(Data data) : calc_data(data) {};
	virtual ~Solver() = default;

	Data calc_data;

private:
	Solver() = delete;

	Eigen::SparseMatrix <double> K;
	Eigen::SparseVector <double> F;
	// Eigen::SparseMatrix <double> M; dyyn
	// Eigen::SparseMatrix <double> C; dyn

	Eigen::VectorX <double> U; // displacement

};