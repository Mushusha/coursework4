#include "Node.h"

void Node::set_constraints(int comp, double value) {
	constraints.insert({ comp, value });
}
