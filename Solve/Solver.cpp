#include "Solver.h"

Solver::Solver(Data& data) : calc_data(data) {
	fillGlobalK();
	fillGlobalF();
	fillConstraints();
}
Solver::Solver(const Solver& other)
	: calc_data(other.calc_data),
	K(other.K),
	F(other.F),
	U(other.U) {
}
Solver& Solver::operator=(const Solver& other) {
	if (this != &other) {
		calc_data = other.calc_data;
		K = other.K;
		F = other.F;
		U = other.U;
	}
	return *this;
}
Solver::Solver(Solver&& other) noexcept
	: calc_data(std::move(other.calc_data)),
	K(std::move(other.K)),
	F(std::move(other.F)),
	U(std::move(other.U)) {
}
Solver& Solver::operator=(Solver&& other) noexcept {
	if (this != &other) {
		calc_data = std::move(other.calc_data);
		K = std::move(other.K);
		F = std::move(other.F);
		U = std::move(other.U);
	}
	return *this;
}

void Solver::solve() {
	logger& log = logger::log();
	log.print("Start solving");

	calcDisp();
	dispToElem();
	dispToNode();
	calcStrain();
	calcStress();

	log.print("Solve done");
}

void Solver::fillGlobalK() {
	logger& log = logger::log();
	log.print("Start filling stiffness matrix");
	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();
	int dim = calc_data.dim;

	K.resize(dim * nodes_count, dim * nodes_count);
	
	std::vector<Eigen::Triplet<std::complex<double>>> tripl_vec;
	
	const int PARALLEL_THRESHOLD = 100;
	
	if (elems_count < PARALLEL_THRESHOLD) {
		int avg_nodes_per_elem = elems_count > 0 ? calc_data.get_elem(0)->nodes_count() : 8;
		tripl_vec.reserve(elems_count * avg_nodes_per_elem * avg_nodes_per_elem * dim * dim);
		
		for (int i = 0; i < elems_count; i++) {
			auto elem = calc_data.get_elem(i);
			Eigen::MatrixXcd loc_k = elem->localK();
			int elem_nodes = elem->nodes_count();
			const auto& elem_nodes_vec = elem->get_node();
			
			for (int j = 0; j < elem_nodes * dim; j++)
				for (int k = 0; k < elem_nodes * dim; k++) {
					tripl_vec.push_back(Eigen::Triplet<std::complex<double>>(
						dim * (elem_nodes_vec[j / dim] - 1) + j % dim,
						dim * (elem_nodes_vec[k / dim] - 1) + k % dim, 
						loc_k(j, k)));
				}
		}
	}
	else {
		unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
		std::vector<std::vector<Eigen::Triplet<std::complex<double>>>> chunk_triplets(num_threads);
		
		int chunk_size = (elems_count + num_threads - 1) / num_threads;
		
		std::vector<unsigned int> chunk_indices(num_threads);
		std::iota(chunk_indices.begin(), chunk_indices.end(), 0);
		
		int avg_nodes_per_elem = calc_data.get_elem(0)->nodes_count();
		int triplets_per_elem = avg_nodes_per_elem * avg_nodes_per_elem * dim * dim;
		
		std::for_each(std::execution::par, chunk_indices.begin(), chunk_indices.end(), [&](unsigned int chunk_id) {
			int start = chunk_id * chunk_size;
			int end = std::min(start + chunk_size, elems_count);
			
			if (start >= end) return;
			
			auto& local_triplets = chunk_triplets[chunk_id];
			local_triplets.reserve((end - start) * triplets_per_elem);
			
		for (int i = start; i < end; i++) {
			auto elem = calc_data.get_elem(i);
			Eigen::MatrixXcd loc_k = elem->localK();
			int elem_nodes = elem->nodes_count();
			const auto& elem_nodes_vec = elem->get_node();
			
			for (int j = 0; j < elem_nodes * dim; j++)
				for (int k = 0; k < elem_nodes * dim; k++) {
					local_triplets.push_back(Eigen::Triplet<std::complex<double>>(
						dim * (elem_nodes_vec[j / dim] - 1) + j % dim,
						dim * (elem_nodes_vec[k / dim] - 1) + k % dim, 
						loc_k(j, k)));
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

	K.setFromTriplets(tripl_vec.begin(), tripl_vec.end());
	log.print("End filling stiffness matrix");
}

void Solver::fillGlobalF(double mult) {
	// logger& log = logger::log();
	// log.print("Start filling right vector");

	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();
	int dim = calc_data.dim;

	if (F.size() != dim * nodes_count) {
		F.resize(dim * nodes_count);
	}
	F.setZero();

	const int PARALLEL_THRESHOLD = 100;
	if (elems_count >= PARALLEL_THRESHOLD) {
		unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
		std::vector<Eigen::VectorXcd> thread_local_F(num_threads, Eigen::VectorXcd::Zero(dim * nodes_count));
		
		std::vector<unsigned int> chunk_indices(num_threads);
		std::iota(chunk_indices.begin(), chunk_indices.end(), 0);
		int chunk_size = (elems_count + num_threads - 1) / num_threads;
		
		std::for_each(std::execution::par, chunk_indices.begin(), chunk_indices.end(), [&](unsigned int chunk_id) {
			int start = chunk_id * chunk_size;
			int end = std::min(start + chunk_size, elems_count);
			
			if (start >= end) return;
			
			auto& local_F = thread_local_F[chunk_id];
			
			for (int i = start; i < end; i++) {
				auto elem = calc_data.get_elem(i);
				std::vector<double> loc_f = elem->localF(mult);
				if (loc_f.size() != 0) {
					int elem_nodes = elem->nodes_count();
					const auto& elem_nodes_vec = elem->get_node();
					for (int j = 0; j < elem_nodes * dim; j++) {
						int node_idx = elem_nodes_vec[j / dim] - 1;
						local_F(dim * node_idx + j % dim) += loc_f[j];
					}
				}
			}
		});
		
		// Merge thread-local results
		for (auto& local_F : thread_local_F) {
			F += local_F;
		}
	}
	else {
		for (int i = 0; i < elems_count; i++) {
			auto elem = calc_data.get_elem(i);
			std::vector<double> loc_f = elem->localF(mult);
			if (loc_f.size() != 0) {
				int elem_nodes = elem->nodes_count();
				const auto& elem_nodes_vec = elem->get_node();
				for (int j = 0; j < elem_nodes * dim; j++) {
					int node_idx = elem_nodes_vec[j / dim] - 1;
					F(dim * node_idx + j % dim) += loc_f[j];
				}
			}
		}
	}

	std::vector<size_t> ind(nodes_count);
	std::iota(ind.begin(), ind.end(), 0);

	std::for_each(std::execution::par, ind.begin(), ind.end(), [&](size_t i) {
		for (auto pair : calc_data.get_node(i)->load) {
			F(dim * i + pair.first) += mult * pair.second;
		}
		});
	// log.print("End filling right vector");
}

void Solver::fillConstraints() {
	logger& log = logger::log();
	log.print("Start filling constraints");

	int dim = calc_data.dim;
	int nodes_count = calc_data.nodes_count();
	
	std::vector<std::pair<int, double>> constrained_dofs;
	constrained_dofs.reserve(nodes_count * dim);
	
	std::unordered_set<int> constrained_dof_set;
	
	for (int node = 0; node < nodes_count; node++) {
		for (auto const& c : calc_data.get_node(node)->constraints) {
			int dof = node * dim + c.first;
			constrained_dofs.emplace_back(dof, c.second);
			constrained_dof_set.insert(dof);
		}
	}
	
	if (constrained_dofs.empty()) {
		log.print("End filling constraints (no constraints found)");
		return;
	}
	
	std::vector<int> col_indices(K.outerSize());
	std::iota(col_indices.begin(), col_indices.end(), 0);
	
	std::for_each(std::execution::par, col_indices.begin(), col_indices.end(), [&](int col) {
		for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(K, col); it; ++it) {
			int row = it.row();
			int col_idx = it.col();
			
			if (row != col_idx && 
				(constrained_dof_set.count(row) > 0 || constrained_dof_set.count(col_idx) > 0)) {
				it.valueRef() = 0.0;
			}
		}
	});
	
	std::for_each(std::execution::par, constrained_dofs.begin(), constrained_dofs.end(), [&](const auto& dof_pair) {
		int dof = dof_pair.first;
		double constraint_value = dof_pair.second;
		
		std::complex<double> diagonal_value = K.coeff(dof, dof);
		
		F(dof) = constraint_value * diagonal_value;
	});
	
	log.print("End filling constraints");
}

void Solver::addToGlobalK(int first_index, int second_index, double value) {
	Eigen::Triplet <double> tripl(first_index, second_index, value);
	K.setFromTriplets(&tripl, &tripl + 1);
}

void Solver::addToGlobalF(int index, double value) {
	F(index) = value;
}

void Solver::zeroDiagonalCheck() {
	int total_dofs = calc_data.nodes_count() * calc_data.dim;
	int dim = calc_data.dim;
	
	std::vector<size_t> dof_indices(total_dofs);
	std::iota(dof_indices.begin(), dof_indices.end(), 0);
	
	std::for_each(std::execution::par, dof_indices.begin(), dof_indices.end(), [&](size_t i) {
		if (K.coeffRef(i, i) == std::complex<double>(0, 0)) {
			throw runtime_error("Zero on diagonal: node " + std::to_string(static_cast<int>(i / dim + 1)) + " dof " + std::to_string(static_cast<int>(i % dim)));
		}
	});
}

void Solver::dispToElem() {
	int dim = calc_data.dim;
	int elems_count = calc_data.elements_count();

	std::vector<size_t> elem_indices(elems_count);
	std::iota(elem_indices.begin(), elem_indices.end(), 0);
	
	std::for_each(std::execution::par, elem_indices.begin(), elem_indices.end(), [&](size_t elem) {
		auto el = calc_data.get_elem(elem);
		int elem_nodes = el->nodes_count();
		el->results.resize(elem_nodes);
		el->displacements.resize(dim * elem_nodes);
		const auto& elem_nodes_vec = el->get_node();
		
		for (int node = 0; node < elem_nodes; node++) {
			int global_node_idx = elem_nodes_vec[node] - 1;
			int base_idx = dim * global_node_idx;
			el->displacements[dim * node] = U(base_idx).real();
			el->displacements[dim * node + 1] = U(base_idx + 1).real();
			if (dim == 3)
				el->displacements[dim * node + 2] = U(base_idx + 2).real();
		}
	});
}

void Solver::dispToNode() {
	int nodes_count = calc_data.nodes_count();
	int dim = calc_data.dim;
	
	std::vector<size_t> node_indices(nodes_count);
	std::iota(node_indices.begin(), node_indices.end(), 0);
	
	std::for_each(std::execution::par, node_indices.begin(), node_indices.end(), [&](size_t node) {
		calc_data.get_node(node)->results[DISPLACEMENT].clear();
		for (int i = 0; i < dim; i++)
			calc_data.get_node(node)->set_result(U(dim * node + i).real(), DISPLACEMENT);
	});
}

void Solver::calcStrain() {
	int elems_count = calc_data.elements_count();
	std::vector<size_t> elem_indices(elems_count);
	std::iota(elem_indices.begin(), elem_indices.end(), 0);
	
	std::for_each(std::execution::par, elem_indices.begin(), elem_indices.end(), [&](size_t elem) {
		auto el = calc_data.get_elem(elem);
		for (int node = 0; node < el->nodes_count(); node++) {
			el->results[node][STRAIN] = el->B(el->gaussPoint(KSI, node), el->gaussPoint(ETA, node),
				el->gaussPoint(ZETA, node)) * el->displacements;
		}
	});
	
	logger& log = logger::log();
	log.print("Calculate Strain");
}

void Solver::calcStress() {
	int elems_count = calc_data.elements_count();
	int dim = calc_data.dim;
	std::vector<size_t> elem_indices(elems_count);
	std::iota(elem_indices.begin(), elem_indices.end(), 0);
	
	std::for_each(std::execution::par, elem_indices.begin(), elem_indices.end(), [&](size_t elem) {
		auto el = calc_data.get_elem(elem);
		for (int node = 0; node < el->nodes_count(); node++) {
			el->results[node][STRESS] = el->D * el->results[node][STRAIN];
			if (dim == 2)
				el->results[node][STRAIN][Comp2D::XY] /= 2;
		}
	});
	
	logger& log = logger::log();
	log.print("Calculate Stress");
}
