#pragma once

#include "Parser.h"
#include "Data.h"
#include "Solver.h"
#include "Smoothing/Smoothing.h"
#include "Smoothing/Interpolation.h"
#include "Output/Output.h"
#include "Output/VTUWriter.h"
#include "Fabric.h"

class Calculate {
public:
	Calculate(std::string name) : file(name) {}

	void Solve();

private:
	Calculate();

	std::string file;
};