#pragma once

#include <map>
#include <vector>
#include <array>
#include <utility>

#include "Enums.h"

class Node {
private:
	int id;
	double x;
	double y;
	double z;
	std::array<double, 3> displacement;
	std::array<double, 6> strain; // agreed
	std::array<double, 6> stress; // agreed

public:
	Node() {};
	Node(int id, std::array <double, 3> coords) {
		this->id = id;
		this->x = coords[0];
		this->y = coords[1];
		this->z = coords[2];
	};
	virtual ~Node() = default;

	std::map <int, double> constraints; // first - component, second - value
	void set_constraints(int comp, double value); // need ??
	void set_displacement(double u, int i);
	void set_strain(double epsilon, int i);
	void set_stress(double sigma, int i);

	int getID() { return id; }
	double getX() { return x; }
	double getY() { return y; }
	double getZ() { return z; }
	double get_displacement(int i) { return displacement[i]; }
	double get_strain(int i) { return strain[i]; }
	double get_stress(int i) { return stress[i]; }
};