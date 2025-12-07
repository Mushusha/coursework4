#include "Smoothing.h"

Smoothing::Smoothing(Data& data, std::vector<ResType> type)
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
	
	std::vector<Eigen::Triplet<double>> tripl_vec;
	
	const int PARALLEL_THRESHOLD = 100;
	
	if (elems_count < PARALLEL_THRESHOLD) {
		for (int i = 0; i < elems_count; i++) {
			auto elem = calc_data.get_elem(i);
			Eigen::MatrixXd loc_c = elem->localC();
			int elem_nodes = elem->nodes_count();
			
			for (int j = 0; j < elem_nodes; j++)
				for (int k = 0; k < elem_nodes; k++) {
					tripl_vec.push_back(Eigen::Triplet<double>(
						elem->get_node(j) - 1, 
						elem->get_node(k) - 1, 
						loc_c(j, k)));
				}
		}
	}
	else {
		unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
		std::vector<std::vector<Eigen::Triplet<double>>> chunk_triplets(num_threads);
		
		int chunk_size = (elems_count + num_threads - 1) / num_threads;
		
		std::vector<unsigned int> chunk_indices(num_threads);
		std::iota(chunk_indices.begin(), chunk_indices.end(), 0);
		
		std::for_each(std::execution::par, chunk_indices.begin(), chunk_indices.end(), [&](unsigned int chunk_id) {
			int start = chunk_id * chunk_size;
			int end = std::min(start + chunk_size, elems_count);
			
			if (start >= end) return;
			
			auto& local_triplets = chunk_triplets[chunk_id];
			
			for (int i = start; i < end; i++) {
				auto elem = calc_data.get_elem(i);
				Eigen::MatrixXd loc_c = elem->localC();
				int elem_nodes = elem->nodes_count();
				
				for (int j = 0; j < elem_nodes; j++)
					for (int k = 0; k < elem_nodes; k++) {
						local_triplets.push_back(Eigen::Triplet<double>(
							elem->get_node(j) - 1, 
							elem->get_node(k) - 1, 
							loc_c(j, k)));
					}
			}
		});
		
		size_t total_size = 0;
		for (const auto& tv : chunk_triplets) total_size += tv.size();
		tripl_vec.reserve(total_size);
		
		for (auto& tv : chunk_triplets) {
			tripl_vec.insert(tripl_vec.end(), 
				std::make_move_iterator(tv.begin()), 
				std::make_move_iterator(tv.end()));
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
			value.push_back(calc_data.get_elem(i)->results[j][type](comp).real());
		std::vector<double> loc_r = calc_data.get_elem(i)->localR(value);
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count(); j++)
			R.coeffRef(calc_data.get_elem(i)->get_node(j) - 1) += loc_r[j];
	}
}

void Smoothing::solve() {
	logger& log = logger::log();
	log.print("Start smoothing");
	for (auto t : type) {
		for (int comp = 0; comp < numComp(t, calc_data.dim); comp++) {
			fillGlobalR(t, comp);

			Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
			solver.compute(C);

			if (solver.info() != Eigen::Success)
				throw runtime_error("Error in C");

			Eigen::MatrixXd Result;
			Result = solver.solve(R);

			for (int i = 0; i < calc_data.nodes_count(); i++) {
				calc_data.get_node(i)->set_result(Result(i), t);
			}
		}
	}
	log.print("End smoothing");
}
