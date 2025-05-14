#include "Dynamics.h"

Dynamics::Dynamics(Data& data) : Solver(data) {
	fillGlobalM();
	updateM();
	// ...
	calcDelta_t(data);

	filename = data.filename;

	iter_count = data.max_time / delta_t; //std::min(static_cast<int>(data.max_time / delta_t), data.max_iter);
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
	std::vector <Eigen::Triplet <std::complex<double>>> tripl_vec;
	for (int i = 0; i < elems_count; i++) {
		Eigen::MatrixXcd loc_m = calc_data.get_elem(i)->localM();
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count() * dim; j++)
			for (int k = 0; k < calc_data.get_elem(i)->nodes_count() * dim; k++) {
				Eigen::Triplet <std::complex<double>> trpl(dim * (calc_data.get_elem(i)->get_node(j / dim) - 1) + j % dim, dim * (calc_data.get_elem(i)->get_node(k / dim) - 1) + k % dim, loc_m(j, k));
				tripl_vec.push_back(trpl);
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

	for (int row = 0; row < M.outerSize(); ++row) {
		std::complex<double> rowSum = 0;
		for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(M, row); it; ++it)
			rowSum += it.value();

		if (row < M.cols())
			result.insert(row, row) = rowSum;
	}

	result.makeCompressed();

	M = result;
	log.print("End updating mass matrix");
}

void Dynamics::calcDelta_t(Data& data) {
	double h_min = data.get_elem(0)->Volume();
	double v_max = 0;
	for (int i = 0; i < data.elements_count() - data.num_inf_elems; i++) {
		h_min = (h_min > data.get_elem(i)->Volume()) ? data.get_elem(i)->Volume() : h_min;
		double nu = data.get_elem(i)->get_nu();
		double E = data.get_elem(i)->get_E();
		double v_p = std::sqrt((E * nu / ((1 + nu) * (1 - 2 * nu)) + E / (2 * (1 + nu))) / data.get_elem(i)->get_rho());
		v_max = (v_p > v_max) ? v_p : v_max;
	}
	delta_t = 0.8 * h_min / v_max;
	delta_t = 0.00099157 / 2;
}

void Dynamics::U_curr(Eigen::VectorXcd U_prev, Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev) {
	U = U_prev + V_prev * delta_t + A_prev * pow(delta_t, 2) / 2;
}

void Dynamics::V_curr(Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev) {
	V = V_prev + A_prev * delta_t * (1 - beta1) + A * delta_t * beta1;
}

void Dynamics::A_curr(Eigen::VectorXcd U_prev, Eigen::VectorXcd V_prev, Eigen::VectorXcd A_prev) {
	Eigen::SparseMatrix <std::complex<double>> M1 = (1 + alpha * delta_t * beta1) * M;
	Eigen::VectorXcd F1 = F - alpha * M *
		(U_prev + A_prev * delta_t * (1 - beta1)) -
		K * (U_prev + V_prev * delta_t + A_prev * pow(delta_t, 2) / 2);

	A.resize(A_0.size());
	A.setZero();
	//for (int i = 0; i < A.size(); i++)
	//	A(i) = F1(i) / M1.coeffRef(i, i);

	std::vector<int> ind(A.size());
	std::iota(ind.begin(), ind.end(), 0);
	
	std::for_each(std::execution::par, ind.begin(), ind.end(), [&](int i) {
		A(i) = F1(i) / M1.coeffRef(i, i);
	});
}

void Dynamics::calcDisp() {
	Eigen::VectorXcd U_prev = U_0;
	Eigen::VectorXcd V_prev = V_0;
	Eigen::VectorXcd A_prev = A_0;

	for (int i = 0; i < iter_count; i++) {
		fillGlobalF(berlage(omega, Amp, delta_t *  i));
		
		A_curr(U_prev, V_prev, A_prev);
		V_curr(V_prev, A_prev);
		U_curr(U_prev, V_prev, A_prev);

		A_prev = A;
		V_prev = V;
		U_prev = U;

		if (i % iter_res_output == 0) {
			dispToNode();

			VTUWriter vtu(calc_data.get_elements(), calc_data.get_nodes()); // change
			vtu.write(filename + "_" + std::to_string(i / iter_res_output) + ".vtu"); // in directory
		}
	}
	A_curr(U_prev, V_prev, A_prev);
	V_curr(V_prev, A_prev);
	U_curr(U_prev, V_prev, A_prev);
}
