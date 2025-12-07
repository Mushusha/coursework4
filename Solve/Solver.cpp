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
		for (int i = 0; i < elems_count; i++) {
			auto elem = calc_data.get_elem(i);
			Eigen::MatrixXcd loc_k = elem->localK();
			int elem_nodes = elem->nodes_count();
			
			for (int j = 0; j < elem_nodes * dim; j++)
				for (int k = 0; k < elem_nodes * dim; k++) {
					tripl_vec.push_back(Eigen::Triplet<std::complex<double>>(
						dim * (elem->get_node(j / dim) - 1) + j % dim,
						dim * (elem->get_node(k / dim) - 1) + k % dim, 
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
				
				for (int j = 0; j < elem_nodes * dim; j++)
					for (int k = 0; k < elem_nodes * dim; k++) {
						local_triplets.push_back(Eigen::Triplet<std::complex<double>>(
							dim * (elem->get_node(j / dim) - 1) + j % dim,
							dim * (elem->get_node(k / dim) - 1) + k % dim, 
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
	logger& log = logger::log();
	log.print("Start filling right vector");

	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();
	int dim = calc_data.dim;

	F.resize(dim * nodes_count);

	for (int i = 0; i < elems_count; i++) {
		std::vector<double> loc_f = calc_data.get_elem(i)->localF(mult);
		if (loc_f.size() != 0)
			for (int j = 0; j < calc_data.get_elem(i)->nodes_count() * dim; j++)
				F.coeffRef(dim * (calc_data.get_elem(i)->get_node(j / dim) - 1) + j % dim) += loc_f[j];
	}

	//for (int i = 0; i < calc_data.nodes_count(); i++)
	//	for (auto pair : calc_data.get_node(i)->load)
	//		F.coeffRef(dim * i + pair.first) += mult * pair.second;

	std::vector<size_t> ind(calc_data.nodes_count());
	std::iota(ind.begin(), ind.end(), 0);

	std::for_each(std::execution::par, ind.begin(), ind.end(), [&](size_t i) {
		for (auto pair : calc_data.get_node(i)->load) {
			F.coeffRef(dim * i + pair.first) += mult * pair.second;
		}
		});
	log.print("End filling right vector");
}

void Solver::fillConstraints() {
	logger& log = logger::log();
	log.print("Start filling constraints");

	int dim = calc_data.dim;

	// c.first - comp, c.second - value 
	for (int node = 0; node < calc_data.nodes_count(); node++)
		for (auto const& c : calc_data.get_node(node)->constraints)
			for (int i = 0; i < K.outerSize(); i++) {
				for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(K, i); it; ++it)
					if (((it.row() == node * dim + c.first) ||
						 (it.col() == node * dim + c.first)) && (it.row() != it.col())) {
						it.valueRef() = 0.0;
					}
					else if (((it.row() == node * dim + c.first) ||
							  (it.col() == node * dim + c.first)) && (it.row() == it.col()))
						F.coeffRef(it.row()) = c.second * it.value();
			}
	log.print("End filling constraints");
}

void Solver::addToGlobalK(int first_index, int second_index, double value) {
	Eigen::Triplet <double> tripl(first_index, second_index, value);
	K.setFromTriplets(&tripl, &tripl + 1);
}

void Solver::addToGlobalF(int index, double value) {
	F.coeffRef(index) = value;
}

void Solver::zeroDiagonalCheck() {
	for (int i = 0; i < calc_data.nodes_count() * calc_data.dim; i++)
		if (K.coeffRef(i, i) == std::complex<double>(0, 0)) {
			throw runtime_error("Zero on diagonal: node " + std::to_string(static_cast<int>(i / calc_data.dim + 1)) + " dof " + std::to_string(static_cast<int>(i % calc_data.dim)));
			break;
		}
}

void Solver::dispToElem() {
	int dim = calc_data.dim;
	int elems_count = calc_data.elements_count();

	std::vector<size_t> elem_indices(elems_count);
	std::iota(elem_indices.begin(), elem_indices.end(), 0);
	
	std::for_each(std::execution::par, elem_indices.begin(), elem_indices.end(), [&](size_t elem) {
		auto el = calc_data.get_elem(elem);
		el->results.resize(el->nodes_count());
		el->displacements.resize(dim * el->nodes_count());
		
		for (int node = 0; node < el->nodes_count(); node++) {
			el->displacements[dim * node] = U(dim * (el->get_node(node) - 1)).real();
			el->displacements[dim * node + 1] = U(dim * (el->get_node(node) - 1) + 1).real();
			if (dim == 3)
				el->displacements[dim * node + 2] = U(dim * (el->get_node(node) - 1) + 2).real();
		}
	});
}

void Solver::dispToNode() {
	for (int node = 0; node < calc_data.nodes_count(); node++) {
		calc_data.get_node(node)->results[DISPLACEMENT].clear();
		for (int i = 0; i < calc_data.dim; i++)
			calc_data.get_node(node)->set_result(U(calc_data.dim * node + i).real(), DISPLACEMENT);
	}
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
