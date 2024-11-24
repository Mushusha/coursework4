#include "Node.h"

void Node::set_constraints(int comp, double value) {
	constraints.insert({ comp, value });
}

void Node::set_displacement(double u, int i) {
	this->displacement[i] = u;
}

void Node::set_strain(double epsilon, int i) {
	this->strain[i] = epsilon;
}

void Node::set_stress(double sigma, int i) {
	this->stress[i] = sigma;
}
