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

double Element::get_coord(int loc_node, int comp) {
	if (comp == 0) return x[loc_node];
	else if (comp == 1) return y[loc_node];
	else if (comp == 2) return z[loc_node];
	else {
		throw runtime_error("error: incorrect comp");
		return -1;
	}
}

void Element::set_load(int type, int apply_to, std::array<double, 6> value) {
	if (type == PRESSURE)
		set_pressure(apply_to, value[0]);
	else
		for (int i = 0; i < value.size(); i++) {
			throw runtime_error("Error load type");
			//std::pair <int, int> pair(apply_to, i);
			//load.insert({ pair, value[i] });
		}
}
// need in elems
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

Eigen::MatrixXd Element::planeStrainD() {

	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 3);
	D(0, 0) = 1;
	D(0, 1) = Poisson / (1 - Poisson);
	D(1, 0) = Poisson / (1 - Poisson);
	D(1, 1) = 1;
	D(2, 2) = (1 - 2 * Poisson) / (2 * (1 - Poisson)) ;

	D *= Young * (1 - Poisson) / ((1 + Poisson) * (1 - 2 * Poisson));
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