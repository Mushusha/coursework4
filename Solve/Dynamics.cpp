#include "Dynamics.h"
#include "MathMV.h"

Dynamics::Dynamics(Data& data) : Solver(data) {
	fillGlobalM();
	updateM();
	fillGlobalDamping();
	updateDamping();
	// ...
	calcDelta_t(data);

	filename = data.filename;

	iter_count = std::min(static_cast<int>(data.max_time / delta_t), data.max_iter);
	iter_res_output = data.iter_res_output;
	beta1 = 0.5;
	alpha = data.damping;

	Amp = data.Amp;
	omega = data.omega;

	U_0.resize(data.dim * data.nodes_count());
	U_0.setZero();
	V_0.resize(data.dim * data.nodes_count());
	V_0.setZero();
	A_0.resize(data.dim * data.nodes_count());
	A_0.setZero();
}

void Dynamics::fillGlobalM() {
	logger& log = logger::log();
	log.print("Start filling mass matrix");
	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();
	int dim = calc_data.dim;

	M.resize(dim * nodes_count, dim * nodes_count);
	
	std::vector<Eigen::Triplet<std::complex<double>>> tripl_vec;
	
	const int PARALLEL_THRESHOLD = 100;
	
	if (elems_count < PARALLEL_THRESHOLD) {
		int avg_nodes_per_elem = elems_count > 0 ? calc_data.get_elem(0)->nodes_count() : 8;
		tripl_vec.reserve(elems_count * avg_nodes_per_elem * avg_nodes_per_elem * dim * dim);
		
		for (int i = 0; i < elems_count; i++) {
			auto elem = calc_data.get_elem(i);
			Eigen::MatrixXcd loc_m = elem->localM();
			int elem_nodes = elem->nodes_count();
			const auto& elem_nodes_vec = elem->get_node();
			
			for (int j = 0; j < elem_nodes * dim; j++)
				for (int k = 0; k < elem_nodes * dim; k++) {
					tripl_vec.push_back(Eigen::Triplet<std::complex<double>>(
						dim * (elem_nodes_vec[j / dim] - 1) + j % dim,
						dim * (elem_nodes_vec[k / dim] - 1) + k % dim, 
						loc_m(j, k)));
				}
		}
	}
	else {
		unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
		std::vector<std::vector<Eigen::Triplet<std::complex<double>>>> chunk_triplets(num_threads);
		
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
				Eigen::MatrixXcd loc_m = elem->localM();
				int elem_nodes = elem->nodes_count();
				const auto& elem_nodes_vec = elem->get_node();
				
				for (int j = 0; j < elem_nodes * dim; j++)
					for (int k = 0; k < elem_nodes * dim; k++) {
						local_triplets.push_back(Eigen::Triplet<std::complex<double>>(
							dim * (elem_nodes_vec[j / dim] - 1) + j % dim,
							dim * (elem_nodes_vec[k / dim] - 1) + k % dim, 
							loc_m(j, k)));
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

	M.setFromTriplets(tripl_vec.begin(), tripl_vec.end());
	log.print("End filling mass matrix");
}

void Dynamics::updateM() {
	logger& log = logger::log();
	log.print("Start updating mass matrix");

	Eigen::SparseMatrix<std::complex<double>> result(M.rows(), M.cols());
	result.reserve(Eigen::VectorXi::Constant(M.rows(), 1));

	std::vector<int> row_indices(M.outerSize());
	std::iota(row_indices.begin(), row_indices.end(), 0);
	
	std::for_each(std::execution::par, row_indices.begin(), row_indices.end(), [&](int row) {
		std::complex<double> rowSum = 0;
		for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(M, row); it; ++it)
			rowSum += it.value();

		if (row < M.cols())
			result.insert(row, row) = rowSum;
	});

	result.makeCompressed();

	M = result;
	
	M_diagonal = M.diagonal();
	
	log.print("End updating mass matrix");
}

void Dynamics::fillGlobalDamping() {
	logger& log = logger::log();
	log.print("Start filling global damping matrix");
	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();
	int dim = calc_data.dim;

	Damp.resize(dim * nodes_count, dim * nodes_count);

	std::vector<Eigen::Triplet<std::complex<double>>> tripl_vec;
	tripl_vec.reserve(static_cast<size_t>(elems_count) * 64);

	int inf_elem_count = 0;
	double max_damping_value = 0.0;
	double total_damping_contrib = 0.0;
	int non_zero_contributions = 0;

	const int PARALLEL_THRESHOLD = 100;
	if (elems_count >= PARALLEL_THRESHOLD) {
		unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
		std::vector<std::vector<Eigen::Triplet<std::complex<double>>>> chunk_triplets(num_threads);
		std::vector<int> chunk_inf_elem_count(num_threads, 0);
		std::vector<double> chunk_max_damping(num_threads, 0.0);
		std::vector<double> chunk_total_damping(num_threads, 0.0);
		std::vector<int> chunk_non_zero(num_threads, 0);
		
		std::vector<unsigned int> chunk_indices(num_threads);
		std::iota(chunk_indices.begin(), chunk_indices.end(), 0);
		int chunk_size = (elems_count + num_threads - 1) / num_threads;
		
		std::for_each(std::execution::par, chunk_indices.begin(), chunk_indices.end(), [&](unsigned int chunk_id) {
			int start = chunk_id * chunk_size;
			int end = std::min(start + chunk_size, elems_count);
			
			if (start >= end) return;
			
			auto& local_triplets = chunk_triplets[chunk_id];
			auto& local_inf_count = chunk_inf_elem_count[chunk_id];
			auto& local_max = chunk_max_damping[chunk_id];
			auto& local_total = chunk_total_damping[chunk_id];
			auto& local_non_zero = chunk_non_zero[chunk_id];
			
			for (int i = start; i < end; i++) {
				auto elem = calc_data.get_elem(i);
				Eigen::MatrixXd loc_c = elem->localDamping();
				int elem_nodes = elem->nodes_count();

				if (loc_c.rows() != elem_nodes || loc_c.cols() != elem_nodes) continue;

				ElemType elem_type = elem->get_type();
				bool is_inf_elem = (elem_type == INFQUAD || 
				                    elem_type == INFHEX ||
				                    elem_type == INFQUADSEM ||
				                    elem_type == INFHEXSEM);
				
				if (is_inf_elem) {
					local_inf_count++;
					double elem_max = loc_c.maxCoeff();
					double elem_sum = loc_c.sum();
					if (elem_max > local_max) local_max = elem_max;
					local_total += std::abs(elem_sum);
				}

				const auto& elem_nodes_vec = elem->get_node();
				for (int a = 0; a < elem_nodes; ++a) {
					for (int b = 0; b < elem_nodes; ++b) {
						double v = loc_c(a, b);
						if (v == 0.0) continue;
						local_non_zero++;
						for (int comp = 0; comp < dim; ++comp) {
							int row = dim * (elem_nodes_vec[a] - 1) + comp;
							int col = dim * (elem_nodes_vec[b] - 1) + comp;
							local_triplets.emplace_back(row, col, std::complex<double>(v, 0.0));
						}
					}
				}
			}
		});
		
		for (size_t i = 0; i < num_threads; i++) {
			inf_elem_count += chunk_inf_elem_count[i];
			if (chunk_max_damping[i] > max_damping_value) max_damping_value = chunk_max_damping[i];
			total_damping_contrib += chunk_total_damping[i];
			non_zero_contributions += chunk_non_zero[i];
			tripl_vec.insert(tripl_vec.end(), 
				std::make_move_iterator(chunk_triplets[i].begin()), 
				std::make_move_iterator(chunk_triplets[i].end()));
		}
	}
	else {
		for (int i = 0; i < elems_count; i++) {
			auto elem = calc_data.get_elem(i);
			Eigen::MatrixXd loc_c = elem->localDamping();
			int elem_nodes = elem->nodes_count();

			if (loc_c.rows() != elem_nodes || loc_c.cols() != elem_nodes) continue;

			ElemType elem_type = elem->get_type();
			bool is_inf_elem = (elem_type == INFQUAD || 
			                    elem_type == INFHEX ||
			                    elem_type == INFQUADSEM ||
			                    elem_type == INFHEXSEM);
			
			if (is_inf_elem) {
				inf_elem_count++;
				double elem_max = loc_c.maxCoeff();
				double elem_sum = loc_c.sum();
				if (elem_max > max_damping_value) max_damping_value = elem_max;
				total_damping_contrib += std::abs(elem_sum);
			}

			const auto& elem_nodes_vec = elem->get_node();
			for (int a = 0; a < elem_nodes; ++a) {
				for (int b = 0; b < elem_nodes; ++b) {
					double v = loc_c(a, b);
					if (v == 0.0) continue;
					non_zero_contributions++;
					for (int comp = 0; comp < dim; ++comp) {
						int row = dim * (elem_nodes_vec[a] - 1) + comp;
						int col = dim * (elem_nodes_vec[b] - 1) + comp;
						tripl_vec.emplace_back(row, col, std::complex<double>(v, 0.0));
					}
				}
			}
		}
	}

	Damp.setFromTriplets(tripl_vec.begin(), tripl_vec.end());
}

void Dynamics::updateDamping() {
	logger& log = logger::log();
	log.print("Start updating damping matrix (lumping)");
	
	Eigen::SparseMatrix<std::complex<double>> result(Damp.rows(), Damp.cols());
	result.reserve(Eigen::VectorXi::Constant(Damp.rows(), 1));

	int non_zero_rows = 0;
	double max_diagonal = 0.0;
	double min_diagonal = std::numeric_limits<double>::max();
	double total_diagonal = 0.0;

	unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
	std::vector<int> row_indices(Damp.outerSize());
	std::iota(row_indices.begin(), row_indices.end(), 0);
	
	std::vector<int> local_non_zero_rows(num_threads, 0);
	std::vector<double> local_max_diagonal(num_threads, 0.0);
	std::vector<double> local_min_diagonal(num_threads, std::numeric_limits<double>::max());
	std::vector<double> local_total_diagonal(num_threads, 0.0);
	
	std::for_each(std::execution::par, row_indices.begin(), row_indices.end(), [&](int row) {
		std::complex<double> rowSum = 0;
		for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(Damp, row); it; ++it)
			rowSum += it.value();

		if (row < Damp.cols()) {
			double diag_val = rowSum.real();
			if (std::abs(diag_val) > 1e-10) {
				result.insert(row, row) = rowSum;
				
				size_t thread_id = row % num_threads;
				local_non_zero_rows[thread_id]++;
				if (diag_val > local_max_diagonal[thread_id]) local_max_diagonal[thread_id] = diag_val;
				if (diag_val < local_min_diagonal[thread_id]) local_min_diagonal[thread_id] = diag_val;
				local_total_diagonal[thread_id] += std::abs(diag_val);
			}
		}
	});
	
	for (size_t i = 0; i < num_threads; i++) {
		non_zero_rows += local_non_zero_rows[i];
		if (local_max_diagonal[i] > max_diagonal) max_diagonal = local_max_diagonal[i];
		if (local_min_diagonal[i] < min_diagonal) min_diagonal = local_min_diagonal[i];
		total_diagonal += local_total_diagonal[i];
	}

	result.makeCompressed();
	Damp = result;
	Damp_diagonal = Damp.diagonal();
}

void Dynamics::calcDelta_t(Data& data) {
	double h_min = min_edge_length(data.get_elem(0));
	double v_max = 0;
	for (int i = 0; i < data.elements_count() - data.num_inf_elems; i++) {
		double elem_min_edge = min_edge_length(data.get_elem(i));
		h_min = (h_min > elem_min_edge) ? elem_min_edge : h_min;
		double nu = data.get_elem(i)->get_nu();
		double E = data.get_elem(i)->get_E();
		double v_p = std::sqrt((E * nu / ((1 + nu) * (1 - 2 * nu)) + E / (2 * (1 + nu))) / data.get_elem(i)->get_rho());
		v_max = (v_p > v_max) ? v_p : v_max;
	}
	delta_t = 0.1 * h_min / v_max;
}

void Dynamics::U_curr(Eigen::VectorXcd U_prev, Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev) {
	U = U_prev + V_prev * delta_t + A_prev * pow(delta_t, 2) / 2;
}

void Dynamics::V_curr(Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev) {
	V = V_prev + A_prev * delta_t * (1 - beta1) + A * delta_t * beta1;
}

void Dynamics::A_curr(Eigen::VectorXcd U_prev, Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev) {
	Eigen::VectorXcd temp_U = U_prev + A_prev * delta_t * (1 - beta1);
	Eigen::VectorXcd temp_K = U_prev + V_prev * delta_t + A_prev * pow(delta_t, 2) / 2;
	Eigen::VectorXcd temp_V = V_prev + A_prev * delta_t * (1 - beta1);
	
	Eigen::VectorXcd F1;
	F1.noalias() = F - Damp * temp_V - alpha * M * temp_U - K * temp_K;

	A.resize(A_0.size());
	
	const std::complex<double> M1_factor = 1.0 + alpha * delta_t * beta1;
	const std::complex<double> Damp_factor = beta1 * delta_t;
	
	std::vector<int> ind(A.size());
	std::iota(ind.begin(), ind.end(), 0);
	
	std::for_each(std::execution::par, ind.begin(), ind.end(), [&](int i) {
		A(i) = F1(i) / (M1_factor * M_diagonal(i) + Damp_factor * Damp_diagonal(i));
	});
}

void Dynamics::calcDisp() {
	std::filesystem::create_directory(filename);
	std::string full_path = (std::filesystem::path(filename) / filename).string();

	Eigen::VectorXcd U_prev = U_0;
	Eigen::VectorXcd V_prev = V_0;
	Eigen::VectorXcd A_prev = A_0;

	for (int i = 0; i < iter_count; i++) {
		fillGlobalF(berlage(omega, Amp, delta_t * i));

		A_curr(U_prev, V_prev, A_prev);
		V_curr(V_prev, A_prev);
		U_curr(U_prev, V_prev, A_prev);

		A_prev = A;
		V_prev = V;
		U_prev = U;


		if (i % iter_res_output == 0) {
			dispToNode();

			VTUWriter vtu(calc_data.get_elements(), calc_data.get_nodes());
			vtu.write(full_path + "_" + std::to_string(i / iter_res_output) + ".vtu");
		}
	}
	A_curr(U_prev, V_prev, A_prev);
	V_curr(V_prev, A_prev);
	U_curr(U_prev, V_prev, A_prev);
}
