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
	std::vector <std::vector <double>> results; // agreed

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
	void set_result(double value, int type, int comp);
	void set_res_size(int type_size, int dim);

	int getID() { return id; }
	double getX() { return x; }
	double getY() { return y; }
	double getZ() { return z; }
	double get_result(int type, int comp) { return results[type][comp]; }
};