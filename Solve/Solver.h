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

	void solve();

protected:
	Data calc_data;

	Eigen::SparseMatrix <std::complex<double>> K;
	Eigen::SparseVector <std::complex<double>> F;

	Eigen::VectorX <std::complex<double>> U; // displacement

	void fillGlobalF(double mult = 1);

	void addToGlobalK(int first_index, int second_index, double value);
	void addToGlobalF(int index, double value);

	void zeroDiagonalCheck();

	virtual void calcDisp() = 0;

private:
	Solver() = delete;

	void fillGlobalK();

	void fillConstraints();

	void dispToElem();
	void calcStrain();
	void calcStress();
};