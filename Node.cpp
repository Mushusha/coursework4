#include "Node.h"

void Node::set_constraints(int comp, double value) {
	constraints.insert({ comp, value });
}

void Node::set_result(double value, int type, int comp) {
	results[type][comp] = value;
}

void Node::set_res_size(int type_size, int dim) {
	results.resize(COUNT);
	results[DISPLACEMENT].resize(dim);
	if (dim == 2) {
		results[STRAIN].resize(3);
		results[STRESS].resize(3);
	}
	else {
		results[STRAIN].resize(6);
		results[STRESS].resize(6);
	}
}


