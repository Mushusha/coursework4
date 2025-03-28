#include "Output.h"

Output::Output(std::vector<std::vector<double>>& points, 
			   std::vector<std::vector<double>>& values, 
			   ResType type, int dim) 
	: points(points), values(values),
	  type(type), dim(dim) {}

void Output::write() {
	std::ofstream file;

	for (int comp = 0; comp < points.size(); comp++) {
		std::string name = filename(type, comp, dim);
		file.open(name + ".txt");

		//double line_len = sqrt(pow((line_end[0] - line_start[0]), 2) + pow((line_end[1] - line_start[1]), 2));
		for (int i = 0; i < points[comp].size(); i++)
			//file << "(" << line_len / points_count * i << ", " << std::fixed << std::setprecision(12) << values[i] << ")" << std::endl;
			file << "(" << points[i][0] << ", " << std::fixed << std::setprecision(12) << values[comp][i] << ")" << std::endl;
			//file << points[i][0] << std::endl;
			//file << std::fixed << std::setprecision(12) << values[i] << std::endl;

		file.close();
	}
}

std::string Output::filename(ResType type, int comp, int dim) {
	std::string name = "";
	switch (type) {
		case STRAIN:
			name += "strain";
			break;
		case STRESS:
			name += "stress";
			break;
		case DISPLACEMENT:
			name += "displacement";
			break;
		case VELOCITY:
			name += "velocity";
			break;
		case ACCELERATION:
			name += "acceleration";
			break;
		case COUNT:
		default:
			break;
	}
	if (dim == 2) {
		using namespace Comp2D;
			switch (static_cast<Comp>(comp)) {
				case(XX):
					name += "_xx";
					break;
				case(YY):
					name += "_yy";
					break;
				case(XY):
					name += "_xy";
					break;
				default:
					break;
			}
	}
	else if (dim == 3) {
		using namespace Comp3D;
		switch (static_cast<Comp>(comp)) {
			case(XX):
				name += "_xx";
				break;
			case(YY):
				name += "_yy";
				break;
			case(ZZ):
				name += "_zz";
				break;
			case(XY):
				name += "_xy";
				break;
			case(XZ):
				name += "_xz";
				break;
			case(YZ):
				name += "_yz";
				break;
			default:
				break;
		}
	}
	return name;
}
