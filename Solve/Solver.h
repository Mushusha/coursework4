#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/Core"
#include "unsupported/Eigen/IterativeSolvers"

#include "Calc/Data.h"
#include "log1.h"

class Solver {
public:
	Solver(Data& data);
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);
	Solver(Solver&& other) noexcept;
	Solver& operator=(Solver&& other) noexcept;
	virtual ~Solver() = default;

	Data calc_data;

	void solve();

protected:
	Eigen::SparseMatrix <double> K;
	Eigen::SparseVector <double> F;

	Eigen::VectorX <double> U; // displacement

	void fillGlobalF(int n = 0);

	void addToGlobalK(int first_index, int second_index, double value);
	void addToGlobalF(int index, double value);

	void zeroDiagonalCheck();

	virtual void calcDisp() = 0;

private:
	Solver() = delete;

	void fillGlobalK();

	void fillConstraints();

	void displacementToElements();
	void displacementToNodes();
	void calcStrain();
	void calcStress();
};