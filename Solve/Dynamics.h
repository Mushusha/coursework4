#pragma once

#include "Solver.h"

class Dynamics : public Solver { // explicit
public:
	Dynamics(Data& data);
	virtual ~Dynamics() = default;

private:
	Dynamics() = delete;

	Eigen::SparseMatrix <double> M; 
	//Eigen::SparseMatrix <double> C;
	Eigen::VectorX <double> V;
	Eigen::VectorX <double> A;

	void fillGlobalM();
	// void fillGlobalC();

	double beta1;
	double delta_t;
	double alpha;
	int iter_count;

	void calcDelta_t(Data& data);

	Eigen::VectorX <double> U_0;
	Eigen::VectorX <double> V_0;
	Eigen::VectorX <double> A_0;

	void U_curr(Eigen::VectorXd U_prev, Eigen::VectorXd V_prev, Eigen::VectorXd A_prev);
	void V_curr(Eigen::VectorXd V_prev, Eigen::VectorXd A_prev);
	void A_curr(Eigen::VectorXd U_prev, Eigen::VectorXd V_prev, Eigen::VectorXd A_prev);

	void calcDisp() final;
};