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
	std::vector <std::vector <double>> results; // agreed

public:
	Node(int id, std::array<double, 3> coords)
        : id(id), 
		x(coords[0]), 
		y(coords[1]), 
		z(coords[2]) {}
	Node(const Node& other)
	: id(other.id), 
		x(other.x), 
		y(other.y), 
		z(other.z), 
		results(other.results), 
		constraints(other.constraints) {}
	Node& operator=(const Node& other) {
		if (this != &other) {
			id = other.id;
			x = other.x;
			y = other.y;
			z = other.z;
			results = other.results;
			constraints = other.constraints;
		}
		return *this;
	}
	Node(Node&& other) noexcept
		: id(std::exchange(other.id, 0)),
		x(std::exchange(other.x, 0.0)),
		y(std::exchange(other.y, 0.0)),
		z(std::exchange(other.z, 0.0)),
		results(std::move(other.results)),
		constraints(std::move(other.constraints)) {}
	Node& operator=(Node&& other) noexcept {
		if (this != &other) {
			id = std::exchange(other.id, 0);
			x = std::exchange(other.x, 0.0);
			y = std::exchange(other.y, 0.0);
			z = std::exchange(other.z, 0.0);
			results = std::move(other.results);
			constraints = std::move(other.constraints);
		}
		return *this;
	}
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

	void printStress();
};