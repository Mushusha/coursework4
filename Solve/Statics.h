#pragma once

#include "Solver.h"

class Statics : public Solver {
public:
	Statics(Data& data) : Solver(calc_data) {};
	virtual ~Statics() = default;

private:
	Statics() = delete;

	void calcDisp() final;
};