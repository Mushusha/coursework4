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

double Element::Jac(double ksi, double eta, double zeta) {
	return J(ksi, eta, zeta).determinant();
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
