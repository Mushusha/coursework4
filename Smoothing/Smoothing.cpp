#include "Smoothing.h"

Smoothing::Smoothing(Data& data, std::vector<double> start, std::vector<double> end, int count)
	: calc_data(data), line_start(start), line_end(end), points_count(count) {
	fillGlobalC();
}
Smoothing::Smoothing(const Smoothing& other)
    : calc_data(other.calc_data),
    line_start(other.line_start),
    line_end(other.line_end),
    points_count(other.points_count),
    C(other.C),
    R(other.R) {
}
Smoothing& Smoothing::operator=(const Smoothing& other) {
    if (this != &other) {
        calc_data = other.calc_data;
        line_start = other.line_start;
        line_end = other.line_end;
        points_count = other.points_count;
        C = other.C;
        R = other.R;
    }
    return *this;
}
Smoothing::Smoothing(Smoothing&& other)
    : calc_data(std::move(other.calc_data)),
    line_start(std::move(other.line_start)),
    line_end(std::move(other.line_end)),
    points_count(std::move(other.points_count)),
    C(std::move(other.C)),
    R(std::move(other.R)) {
}
Smoothing& Smoothing::operator=(Smoothing&& other) {
    if (this != &other) {
        calc_data = std::move(other.calc_data);
        line_start = std::move(other.line_start);
        line_end = std::move(other.line_end);
        points_count = std::move(other.points_count);
        C = std::move(other.C);
        R = std::move(other.R);
    }
    return *this;
}

void Smoothing::fillGlobalC() {
	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();

	C.resize(nodes_count, nodes_count);
	std::vector <Eigen::Triplet <double>> tripl_vec;
	for (int i = 0; i < elems_count; i++) {
		Eigen::MatrixXd loc_c = calc_data.get_elem(i)->localC();
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count(); j++)
			for (int k = 0; k < calc_data.get_elem(i)->nodes_count(); k++) {
				Eigen::Triplet <double> trpl(calc_data.get_elem(i)->get_nodes(j) - 1, calc_data.get_elem(i)->get_nodes(k) - 1, loc_c(j, k));
				tripl_vec.push_back(trpl);
			}
	}

	C.setFromTriplets(tripl_vec.begin(), tripl_vec.end());
}

void Smoothing::fillGlobalR(int type, int comp) {
	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();

	R.resize(nodes_count);
	for (int i = 0; i < elems_count; i++) {
		std::vector<double> value;
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count(); j++)
			value.push_back(calc_data.get_elem(i)->results[j][type](comp));
		std::vector<double> loc_r = calc_data.get_elem(i)->localR(value);
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count(); j++)
			R.coeffRef(calc_data.get_elem(i)->get_nodes(j) - 1) += loc_r[j];
	}
}
