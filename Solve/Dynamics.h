#pragma once

#include "Solver.h"

class Dynamics : public Solver{
public:
	Dynamics(Data& data);
	virtual ~Dynamics() = default;

private:
	Dynamics() = delete;

	Eigen::SparseMatrix <double> M; 
	Eigen::SparseMatrix <double> C;

	void fillGlobalM();
	// void fillGlobalC();

	void calcDisp() final;
};