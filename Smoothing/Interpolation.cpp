#include "Interpolation.h"

Interpolation::Interpolation(Data& data, std::vector<double> start, std::vector<double> end, int count, ResType type)
	: calc_data(data), type(type), count(count) {
	points.resize(count);
	for (int i = 0; i < count; i++) {
		points[i].resize(start.size());
		for (int j = 0; j < start.size(); j++)
			points[i][j] = (end[j] - start[j]) * i / (count - 1) + start[j];
	}
	values.resize(numComp(type, data.dim));
	for (int i = 0; i < numComp(type, data.dim); i++)
		values[i].resize(count);
}

void Interpolation::solve() {
	for (int p = 0; p < count; p++)
		for (int i = 0; i < calc_data.elements_count(); i++) {
			if (calc_data.get_elem(i)->pointInElem(points[p])) {
				//values[p] = elements[i]->results[0][type][comp];
				std::vector<double> Coord = calc_data.get_elem(i)->coordFF(points[p][0], points[p][1]); // dim == 3
				std::vector<double> N = calc_data.get_elem(i)->FF(Coord[0], Coord[1]);
				for (int node = 0; node < calc_data.get_elem(i)->nodes_count(); node++)
					for (int comp = 0; comp < numComp(type, calc_data.dim); comp++)
						if (isTensor(type))
							values[comp][p] += N[node] * calc_data.get_node(calc_data.get_elem(i)->get_nodes(node) - 1)->get_result(type, comp);
						else if (isVector(type))
							values[comp][p] += N[node] * calc_data.get_elem(i)->displacements(node * calc_data.dim + comp);
					//if (type != DISPLACEMENT)
					//	values[p] += N[node] * elements[i]->results[node][type](comp);
					//else
					//	values[p] += N[node] * elements[i]->displacements(node * dim + comp);
				break;
			}
		}
}
