#include "Dynamics.h"

Dynamics::Dynamics(Data& data) : Solver(data) {
	fillGlobalM();
	updateM();
	// ...
	calcDelta_t(data);
	iter_count = std::min(static_cast<int>(data.max_time / delta_t), data.max_iter);
	beta1 = 0.5;
	alpha = data.damping;

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
	std::vector <Eigen::Triplet <double>> tripl_vec;
	for (int i = 0; i < elems_count; i++) {
		Eigen::MatrixXd loc_m = calc_data.get_elem(i)->localM();
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count() * dim; j++)
			for (int k = 0; k < calc_data.get_elem(i)->nodes_count() * dim; k++) {
				Eigen::Triplet <double> trpl(dim * (calc_data.get_elem(i)->get_nodes(j / dim) - 1) + j % dim, dim * (calc_data.get_elem(i)->get_nodes(k / dim) - 1) + k % dim, loc_m(j, k));
				tripl_vec.push_back(trpl);
			}
	}

	M.setFromTriplets(tripl_vec.begin(), tripl_vec.end());
	log.print("End filling mass matrix");
}

void Dynamics::updateM() {
	logger& log = logger::log();
	log.print("Start updating mass matrix");

	Eigen::SparseMatrix<double> result(M.rows(), M.cols());
	result.reserve(Eigen::VectorXi::Constant(M.rows(), 1));

	for (int row = 0; row < M.outerSize(); ++row) {
		double rowSum = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, row); it; ++it)
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
	for (int i = 0; i < data.elements_count(); i++) {
		h_min = (h_min > data.get_elem(i)->Volume()) ? data.get_elem(i)->Volume() : h_min;
		double nu = data.get_elem(i)->get_nu();
		double E = data.get_elem(i)->get_E();
		double v_p = std::sqrt((E * nu / ((1 + nu) * (1 - 2 * nu)) + E / (2 * (1 + nu))) / data.get_elem(i)->get_rho());
		v_max = (v_p > v_max) ? v_p : v_max;
	}
	delta_t = 0.8 * h_min / v_max;
	delta_t = 0.00171221 / 6.5;
}

void Dynamics::U_curr(Eigen::VectorXd U_prev, Eigen::VectorXd V_prev, Eigen::VectorXd A_prev) {
	U = U_prev + V_prev * delta_t + A_prev * pow(delta_t, 2) / 2;
}

void Dynamics::V_curr(Eigen::VectorXd V_prev, Eigen::VectorXd A_prev) {
	V = V_prev + A_prev * delta_t * (1 - beta1) + A * delta_t * beta1;
}

void Dynamics::A_curr(Eigen::VectorXd U_prev, Eigen::VectorXd V_prev, Eigen::VectorXd A_prev) {
	Eigen::SparseMatrix <double> M1 = (1 + alpha * delta_t * beta1) * M;
	Eigen::VectorX <double> F1 = F - alpha * M *
		(U_prev + A_prev * delta_t * (1 - beta1)) -
		K * (U_prev + V_prev * delta_t + A_prev * pow(delta_t, 2) / 2);

	A.resize(A_0.size());
	A.setZero();
	for (int i = 0; i < A.size(); i++)
		A(i) = F1(i) / M1.coeffRef(i, i);
}

void Dynamics::calcDisp() {
	Eigen::VectorX <double> U_prev = U_0;
	Eigen::VectorX <double> V_prev = V_0;
	Eigen::VectorX <double> A_prev = A_0;

	for (int i = 0; i < iter_count; i++) {
		fillGlobalF(i);

		A_curr(U_prev, V_prev, A_prev);
		V_curr(V_prev, A_prev);
		U_curr(U_prev, V_prev, A_prev);
		
		//printVector(A);
		//std::cout << "-----------------------    A   ----------------------" << std::endl;
		//printVector(V);
		//std::cout << "-----------------------    V   ----------------------" << std::endl;
		//printVector(U);
		//std::cout << "-----------------------    U   ----------------------" << std::endl;
		//std::cout << "-----------------------    " << i << "   ----------------------" << std::endl;


		A_prev = A;
		V_prev = V;
		U_prev = U;
	}
	printVector(U);
	A_curr(U_prev, V_prev, A_prev);
	V_curr(V_prev, A_prev);
	U_curr(U_prev, V_prev, A_prev);
}
