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
	U(std::move(other.U)) {}
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
	std::vector <Eigen::Triplet <double>> tripl_vec;
	for (int i = 0; i < elems_count; i++) {
		Eigen::MatrixXd loc_k = calc_data.get_elem(i)->localK();
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count() * dim; j++)
			for (int k = 0; k < calc_data.get_elem(i)->nodes_count() * dim; k++) {
				Eigen::Triplet <double> trpl(dim * (calc_data.get_elem(i)->get_nodes(j / dim) - 1) + j % dim, dim * (calc_data.get_elem(i)->get_nodes(k / dim) - 1) + k % dim, loc_k(j, k));
				tripl_vec.push_back(trpl);
			}
	}

	K.setFromTriplets(tripl_vec.begin(), tripl_vec.end());
	log.print("End filling stiffness matrix");
}

void Solver::fillGlobalF(int n) {
	logger& log = logger::log();
	log.print("Start filling right vector");

	int nodes_count = calc_data.nodes_count();
	int elems_count = calc_data.elements_count();
	int dim = calc_data.dim;

	F.resize(dim * nodes_count);
	for (int i = 0; i < elems_count; i++) {
		std::vector<double> loc_f = calc_data.get_elem(i)->localF();
		for (int j = 0; j < calc_data.get_elem(i)->nodes_count() * dim; j++)
			F.coeffRef(dim * (calc_data.get_elem(i)->get_nodes(j / dim) - 1) + j % dim) += loc_f[j];
	}
	for (int i = 0; i < calc_data.nodes_count(); i++)
		for (auto pair : calc_data.get_node(i).load)
			F.coeffRef(dim * (i - 1) + pair.first) += pair.second;

	log.print("End filling right vector");
}

void Solver::fillConstraints() {
	logger& log = logger::log();
	log.print("Start filling constraints");

	int dim = calc_data.dim;

	// c.first - comp, c.second - value 
	for (int node = 0; node < calc_data.nodes_count(); node++)
		for (auto const& c : calc_data.get_node(node).constraints)
			for (int i = 0; i < K.outerSize(); i++) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
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
		if (K.coeffRef(i, i) == 0) {
			std::cout << "Zero on diagonal: node " + std::to_string(static_cast<int>(i / calc_data.dim + 1)) + " dof " + std::to_string(static_cast<int>(i % calc_data.dim)) << std::endl;
			break;
		}
}

void Solver::dispToElem() {
	int dim = calc_data.dim;

	for (int elem = 0; elem < calc_data.elements_count(); elem++) {
		calc_data.get_elem(elem)->results.resize(calc_data.get_elem(elem)->nodes_count());
		for (int node = 0; node < calc_data.get_elem(elem)->nodes_count(); node++) {
			//elements[elem]->results[node][DISPLACEMENT].resize(dim);

			//elements[elem]->results[node][DISPLACEMENT][X] = U(dim * (elements[elem]->get_nodes(node) - 1));
			//elements[elem]->results[node][DISPLACEMENT][Y] = U(dim * (elements[elem]->get_nodes(node) - 1) + 1);
			//if (dim == 3)
			//	elements[elem]->results[node][DISPLACEMENT][Z] = U(dim * (elements[elem]->get_nodes(node) - 1) + 2);

			calc_data.get_elem(elem)->displacements.resize(dim * calc_data.get_elem(elem)->nodes_count());
			calc_data.get_elem(elem)->displacements[dim * node] = U(dim * (calc_data.get_elem(elem)->get_nodes(node) - 1));
			calc_data.get_elem(elem)->displacements[dim * node + 1] = U(dim * (calc_data.get_elem(elem)->get_nodes(node) - 1) + 1);
			if (dim == 3)
				calc_data.get_elem(elem)->displacements[dim * node + 2] = U(dim * (calc_data.get_elem(elem)->get_nodes(node) - 1) + 2);
		}
	}
}

void Solver::calcStrain() {
	for (int elem = 0; elem < calc_data.elements_count(); elem++) {
		std::vector <double> ksi = { -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926 };
		std::vector <double> eta = { -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926 };
		std::vector <double> zeta = { -0.57735026918926, -0.57735026918926, -0.57735026918926, -0.57735026918926, 0.57735026918926, 0.57735026918926, 0.57735026918926, 0.57735026918926 };
		//std::vector <double> ksi = { -1, 1, 1, -1, -1, 1, 1, -1 };
		//std::vector <double> eta = { -1, -1, 1, 1, -1, -1, 1, 1 };
		//std::vector <double> zeta = { -1, -1, -1, -1, 1, 1, 1, 1 };

		for (int node = 0; node < calc_data.get_elem(elem)->nodes_count(); node++) {
			calc_data.get_elem(elem)->results[node][STRAIN] = 
				calc_data.get_elem(elem)->B(ksi[node], eta[node], zeta[node]) * calc_data.get_elem(elem)->displacements;
		}
	}
	logger& log = logger::log();
	log.print("Calculate Strain");
}

void Solver::calcStress() {
	for (int elem = 0; elem < calc_data.elements_count(); elem++) {
		for (int node = 0; node < calc_data.get_elem(elem)->nodes_count(); node++) {
			calc_data.get_elem(elem)->results[node][STRESS] = calc_data.get_elem(elem)->D * calc_data.get_elem(elem)->results[node][Tensor::STRAIN];
			if (calc_data.dim == 2)
				calc_data.get_elem(elem)->results[node][STRAIN][Comp2D::XY] /= 2;
		}
	}
	logger& log = logger::log();
	log.print("Calculate Stress");
}
