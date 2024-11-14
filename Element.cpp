#include "Element.h"


void Element::set_coords(std::vector <double> x, std::vector <double> y, std::vector <double> z) {
	for (int i = 0; i < x.size(); i++) {
		this->x.push_back(x[i]);
		this->y.push_back(y[i]);
		this->z.push_back(z[i]);
	}
}

void Element::set_constants(double E, double nu) {
	Young = E;
	Poisson = nu;
}

void Element::set_load(int type, int apply_to, std::array<double, 6> value) {
	if (type == LoadType::PRESSURE)
		set_pressure(apply_to, value[0]);
	else
		for (int i = 0; i < value.size(); i++) {
			std::cout << "Error load type" << std::endl;
			//std::pair <int, int> pair(apply_to, i);
			//load.insert({ pair, value[i] });
		}
}

Eigen::MatrixXd Element::planeStressD() {

	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 3);
	D(0, 0) = 1;
	D(0, 1) = Poisson;
	D(1, 0) = Poisson; 
	D(1, 1) = 1;
	D(2, 2) = (1 - Poisson) / 2;

	D *= Young / (1 - pow(Poisson, 2)); 
	return D;
}

Eigen::MatrixXd Element::threeMatrixD() {
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(6, 6);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 3; j++) {
			if (i == j && i < 3)
				D(i, j) = 1;
			if (i != j && i < 3)
				D(i, j) = Poisson / (1 - Poisson);
			if (i >= 3)
				D(i, i) = (1 - 2 * Poisson) / (2 * (1 - Poisson));
		}

	D *= Young * (1 - Poisson) / ((1 + Poisson) * (1 - 2 * Poisson));
	return D;
}


// createInf ()