#include "Element.h"

Element::Element(int id, ElemType type, std::vector <int> nodes) :
	id(id), type(type), nodes(nodes) {
};
Element::Element(const Element& other)
	: id(other.id),
	type(other.type),
	nodes(other.nodes),
	x(other.x),
	y(other.y),
	z(other.z),
	Young(other.Young),
	Poisson(other.Poisson),
	density(other.density),
	load(other.load),
	results(other.results),
	displacements(other.displacements),
	D(other.D) {
}
Element& Element::operator=(const Element& other) {
	if (this != &other) {
		id = other.id;
		type = other.type;
		nodes = other.nodes;
		x = other.x;
		y = other.y;
		z = other.z;
		Young = other.Young;
		Poisson = other.Poisson;
		density = other.density;
		load = other.load;
		results = other.results;
		displacements = other.displacements;
		D = other.D;
	}
	return *this;
}
Element::Element(Element&& other) noexcept
	: id(std::exchange(other.id, 0)),
	type(std::move(other.type)),
	nodes(std::move(other.nodes)),
	x(std::move(other.x)),
	y(std::move(other.y)),
	z(std::move(other.z)),
	Young(std::exchange(other.Young, 0.0)),
	Poisson(std::exchange(other.Poisson, 0.0)),
	density(std::exchange(other.density, 0.0)),
	load(std::move(other.load)),
	results(std::move(other.results)),
	displacements(std::move(other.displacements)),
	D(std::move(other.D)) {
}
Element& Element::operator=(Element&& other) noexcept {
	if (this != &other) {
		id = std::exchange(other.id, 0);
		type = std::move(other.type);
		nodes = std::move(other.nodes);
		x = std::move(other.x);
		y = std::move(other.y);
		z = std::move(other.z);
		Young = std::exchange(other.Young, 0.0);
		Poisson = std::exchange(other.Poisson, 0.0);
		density = std::exchange(other.density, 0.0);
		load = std::move(other.load);
		results = std::move(other.results);
		displacements = std::move(other.displacements);
		D = std::move(other.D);
	}
	return *this;
}

void Element::set_coords(std::vector <double> x, std::vector <double> y, std::vector <double> z) {
	for (int i = 0; i < x.size(); i++) {
		this->x.push_back(x[i]);
		this->y.push_back(y[i]);
		this->z.push_back(z[i]);
	}
}

void Element::set_constants(double E, double nu, double rho) {
	Young = E;
	Poisson = nu;
	density = rho;
}

double Element::get_coord(int loc_node, int comp) const {
	if (comp == 0) return x[loc_node];
	else if (comp == 1) return y[loc_node];
	else if (comp == 2) return z[loc_node];
	else {
		throw runtime_error("error: incorrect comp");
		return -1;
	}
}

std::complex<double> Element::Jac(double ksi, double eta, double zeta) {
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
