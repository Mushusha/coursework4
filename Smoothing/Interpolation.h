#pragma once

#include "Calc/Data.h"

// DISABLED BECAUSE ADDED VTK WRITER

class Interpolation {
public:
	Interpolation() = default;
	Interpolation(Data& data, std::vector<double> start,
				  std::vector<double> end, int count, ResType type);

	void solve();

	std::vector<std::vector<double>> values;
	std::vector<std::vector<double>> points;

private:
	Data& calc_data;
	ResType type;
	int count;
};