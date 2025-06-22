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

public:
	Node(int id, std::array<double, 3> coords);
	Node(const Node& other);
	Node& operator=(const Node& other);
	Node(Node&& other) noexcept;
	Node& operator=(Node&& other) noexcept;
	virtual ~Node() = default;

	std::map <int, double> load;
	std::map <int, double> constraints; // first - component, second - value
	std::map <ResType, std::vector <double>> results; // agreed tensor

	void set_constraints(int comp, double value); // need ??
	void set_result(double value, ResType type);

	int getID() { return id; }
	double getX() const { return x; }
	double getY() const { return y; }
	double getZ() const { return z; }
	double getCoord(int i) const;
};