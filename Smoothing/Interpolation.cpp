#include "Interpolation.h"

Interpolation::Interpolation(Data& data, std::vector<double> start, std::vector<double> end, int count)
	: calc_data(data), line_start(start), line_end(end), points_count(count) {}

void Interpolation::solve() {
	for (int p = 0; p < points_count; p++)
		for (int i = 0; i < calc_data.elements_count(); i++) {
			if (calc_data.get_elem(i)->pointInElem(points[p])) {
				//values[p] = elements[i]->results[0][type][comp];
				std::vector<double> Coord = calc_data.get_elem(i)->coordFF(points[p][0], points[p][1]); // dim == 3
				std::vector<double> N = calc_data.get_elem(i)->FF(Coord[0], Coord[1]);
				for (int node = 0; node < calc_data.get_elem(i)->nodes_count(); node++)
					values[p] += N[node] * calc_data.get_node(calc_data.get_elem(i)->get_nodes(node) - 1).get_result(type, comp);

					//if (type != DISPLACEMENT)
					//	values[p] += N[node] * elements[i]->results[node][type](comp);
					//else
					//	values[p] += N[node] * elements[i]->displacements(node * dim + comp);
				break;
			}
		}
}

