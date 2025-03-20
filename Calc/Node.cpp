#include "Node.h"

Node::Node(int id, std::array<double, 3> coords)
	: id(id),
	x(coords[0]),
	y(coords[1]),
	z(coords[2]) {
}
Node::Node(const Node& other)
	: id(other.id),
	x(other.x),
	y(other.y),
	z(other.z),
	results(other.results),
	constraints(other.constraints) {
}
Node& Node::operator=(const Node& other) {
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
Node::Node(Node&& other) noexcept
	: id(std::exchange(other.id, 0)),
	x(std::exchange(other.x, 0.0)),
	y(std::exchange(other.y, 0.0)),
	z(std::exchange(other.z, 0.0)),
	results(std::move(other.results)),
	constraints(std::move(other.constraints)) {
}
Node& Node::operator=(Node&& other) noexcept {
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

void Node::set_constraints(int comp, double value) {
	constraints.insert({ comp, value });
}

void Node::set_result(double value, int type, int comp) {
	results[type].push_back(value);
}
