#include "Smoothing.h"

Smoothing::Smoothing(Data& data, ResType type)
	: calc_data(data), type(type) {
	fillGlobalC();
}
Smoothing::Smoothing(const Smoothing& other)
	: calc_data(other.calc_data),
	type(other.type),
	C(other.C) {}
Smoothing& Smoothing::operator=(const Smoothing& other) {
	if (this != &other) {
		calc_data = other.calc_data;
		C = other.C;
		type = other.type;
	}
	return *this;
}
Smoothing::Smoothing(Smoothing&& other) noexcept
	: calc_data(std::move(other.calc_data)),
	type(std::move(other.type)),
	C(std::move(other.C)) {}
Smoothing& Smoothing::operator=(Smoothing&& other) noexcept {
	if (this != &other) {
		calc_data = std::move(other.calc_data);
		type = std::move(other.type);
		C = std::move(other.C);
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

void Smoothing::fillGlobalR(ResType type, int comp) {
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

void Smoothing::solve() {
	logger& log = logger::log();
	log.print("Start smoothing");
	
	for (int comp = 0; comp < numComp(type, calc_data.dim); comp++) {
		fillGlobalR(type, comp);

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(C);

		if (solver.info() != Eigen::Success)
			throw runtime_error("Error in C");

		Eigen::MatrixXd Result;
		Result = solver.solve(R);

		for (int i = 0; i < calc_data.nodes_count(); i++) {
			calc_data.get_node(i)->set_result(Result(i), type);
		}
	}

	log.print("End smoothing");
}
