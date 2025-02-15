#include "Node.h"

void Node::set_constraints(int comp, double value) {
	constraints.insert({ comp, value });
}

void Node::set_result(double value, int type, int comp) {
	results[type][comp] = value;
}

void Node::set_res_size(int type_size, int dim) {
	results.resize(COUNT);
	results[DISPLACEMENT].resize(dim);
	if (dim == 2) {
		results[STRAIN].resize(3);
		results[STRESS].resize(3);
	}
	else {
		results[STRAIN].resize(6);
		results[STRESS].resize(6);
	}
}

void Node::printStress() {
	std::ofstream file_coord_x;
	file_coord_x.open("coord_x.txt", std::ios::app);
	std::ofstream file_coord_y;
	file_coord_y.open("coord_y.txt", std::ios::app);
	//std::ofstream file_coord_z;
	//file_coord_z.open("coord_z.txt", std::ios::app);

	std::ofstream file_stress_xx;
	file_stress_xx.open("all_stress_xx.txt", std::ios::app);
	std::ofstream file_stress_yy;
	file_stress_yy.open("all_stress_yy.txt", std::ios::app);
	std::ofstream file_stress_xy;
	file_stress_xy.open("all_stress_xy.txt", std::ios::app);
	//std::ofstream file_stress_xz;
	//file_stress_xz.open("all_stress_xz.txt", std::ios::app);
	//std::ofstream file_stress_yz;
	//file_stress_yz.open("all_stress_yz.txt", std::ios::app);
	//std::ofstream file_stress_zz;
	//file_stress_zz.open("all_stress_zz.txt", std::ios::app);

	file_coord_x << x << std::endl;
	file_coord_y << y << std::endl;
	//file_coord_z << z << std::endl;

	file_stress_xx << results[STRESS][XX_2D] << std::endl;
	file_stress_yy << results[STRESS][YY_2D] << std::endl;
	file_stress_xy << results[STRESS][XY_2D] << std::endl;
}


