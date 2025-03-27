#pragma once

#include "Calc/Data.h"

class Interpolation {
public:
	Interpolation() = default;
	Interpolation(Data& data, std::vector<double> start,
				  std::vector<double> end, int count);

	void solve();

private:
	Data& calc_data;
	std::vector<double> line_start;
	std::vector<double> line_end;
	int points_count;
};