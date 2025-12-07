#pragma once

#include "Solver.h"
#include "Output/VTUWriter.h"

class Dynamics : public Solver { // explicit
public:
	Dynamics(Data& data);
	Dynamics(const Dynamics& other) : Solver(other) {} // add 
	Dynamics& operator=(const Dynamics& other) {
		if (this != &other) {
			Solver::operator=(other);
		}
		return *this;
	}
	Dynamics(Dynamics&& other) noexcept
		: Solver(std::move(other)) {}
	Dynamics& operator=(Dynamics&& other) noexcept {
		if (this != &other) {
			Solver::operator=(std::move(other));
		}
		return *this;
	}
	virtual ~Dynamics() = default;

private:
	Eigen::SparseMatrix <std::complex<double>> M;
	//Eigen::SparseMatrix <double> C;
	Eigen::VectorX <std::complex<double>> V;
	Eigen::VectorX <std::complex<double>> A;

	void fillGlobalM();
	void updateM();
	// void fillGlobalC();

	std::string filename;

	double beta1;
	double delta_t;
	double alpha;
	int iter_count;
	int iter_res_output;

	double omega;
	double Amp;

	void calcDelta_t(Data& data);

	Eigen::VectorX <double> U_0;
	Eigen::VectorX <double> V_0;
	Eigen::VectorX <double> A_0;

	void U_curr(Eigen::VectorXcd U_prev, Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev);
	void V_curr(Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev);
	void A_curr(Eigen::VectorXcd U_prev, Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev);

	void calcDisp() final;
};