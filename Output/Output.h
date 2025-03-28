#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

#include "Enums.h"

class Output {
public:
	Output() = delete;
	Output(std::vector<std::vector<double>>& points, 
		   std::vector<std::vector<double>>& values, 
		   ResType type, int dim);

	void write();
private:
	std::vector<std::vector<double>> points;
	std::vector<std::vector<double>> values;

	ResType type;
	int dim;

	static std::string filename(ResType type, int comp, int dim);
};