#pragma once
#include "Solver.h"
#include "Statics.h"
#include "Dynamics.h"


class FabricSolver {
public:
	static std::shared_ptr<Solver> createSolver(Data& data) {
		if (data.analisys_type == "dynamic")
			return std::make_shared<Dynamics>(data);

		return std::make_shared<Statics>(data);
	}

private:
	FabricSolver() = default;
	~FabricSolver() = default;
};