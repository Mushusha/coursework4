#pragma once

#include <map>
#include <vector>
#include <array>
#include <utility>
#include <fstream>
#include <iostream>

#include "Enums.h"

class Node {
private:
	Node() = default;

	int id;
	double x;
	double y;
	double z;
	std::array <std::vector <double>, numTensor()> results; // [i] - type, [j] - comp // agreed tensor

public:
	Node(int id, std::array<double, 3> coords);
	Node(const Node& other);
	Node& operator=(const Node& other);
	Node(Node&& other) noexcept;
	Node& operator=(Node&& other) noexcept;
	virtual ~Node() = default;

	std::map <int, double> load;
	std::map <int, double> constraints; // first - component, second - value
	void set_constraints(int comp, double value); // need ??
	void set_result(double value, int type, int comp);

	int getID() { return id; }
	double getX() { return x; }
	double getY() { return y; }
	double getZ() { return z; }
	double get_result(int type, int comp) { return results[type][comp]; }
};