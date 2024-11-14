#include "Data.h"

Data::Data(std::shared_ptr <Parser> p) : parser(p) {
	this->dim = 2; // need read from fc
	create_nodes();
	create_elements();
	create_constants();
	create_constraints();
	create_load();
}

void Data::create_nodes() {
	for (int i = 0; i < parser->mesh.nodes_count; i++) {
		std::array <double, 3> coords;
		for (int j = 0; j < 3; j++)
			coords[j] = parser->mesh.nodes_coord[3 * i + j];
		this->nodes.push_back(Node(parser->mesh.node_id[i], coords));
	}
	logger& log = logger::log();
	log.print("Create nodes");
}

void Data::create_elements() {
	int node_tmp = 0;
	for (int i = 0; i < parser->mesh.elems_count; i++) {
		std::vector<int> elem_nodes;
		std::shared_ptr<Element> elem;
		switch (parser->mesh.elem_types[i]) {
		case ElemType::TRI:
			elem_nodes.resize(3);
			for (int j = 0; j < 3; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 3;
			elem = std::make_shared<triElement>(triElement(parser->mesh.elem_id[i], ElemType::TRI, elem_nodes));
			break;
		case ElemType::QUAD:
			elem_nodes.resize(4);
			for (int j = 0; j < 4; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 4;
			elem = std::make_shared<quadElement>(quadElement(parser->mesh.elem_id[i], ElemType::QUAD, elem_nodes));
			break;
		case ElemType::TETRA:
			elem_nodes.resize(4);
			for (int j = 0; j < 4; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 4;
			elem = std::make_shared<tetraElement>(tetraElement(parser->mesh.elem_id[i], ElemType::TETRA, elem_nodes));
			break;
		case ElemType::HEX:
			elem_nodes.resize(8);
			for (int j = 0; j < 8; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 8;
			elem = std::make_shared<hexElement>(hexElement(parser->mesh.elem_id[i], ElemType::HEX, elem_nodes));
			break;
		default:
			std::cout << "Error: incorrect element " + to_string(i) +
				" type: " + to_string(parser->mesh.elem_types[i]) << std::endl;  //need throw
			break;
		}
		std::vector<double> x, y, z;
		for (int j = 0; j < elem_nodes.size(); j++) {
			x.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1)]);
			y.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 1]);
			z.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 2]); // may be Nodes 
		}
		elem->set_coords(x, y, z);
		this->elements.push_back(elem);
	}
	logger& log = logger::log();
	log.print("Create elements");
}

void Data::create_constants() {
	// need block
	for (int i = 0; i < parser->mesh.elems_count; i++)
		elements[i]->set_constants(parser->material[0].constants[0], parser->material[0].constants[1]);
	logger& log = logger::log();
	log.print("Create constants");
}

void Data::create_constraints() {
	for (int id = 0; id < parser->restraints.size(); id++) {
		auto& restraints = parser->restraints[id];
		for (int node = 0; node < restraints.size; node++)
			for (int i = 0; i < 6; i++)
				if (restraints.flag[i])
					nodes[restraints.apply_to[node] - 1].set_constraints(i, restraints.data[i]); // if exist --> throw
	}
	logger& log = logger::log();
	log.print("Create constraints");
}

void Data::create_load() { // type - pressure
	for (int id = 0; id < parser->load.size(); id++) {
		auto& load = parser->load[id];
		if (load.type == LoadType::PRESSURE) 		// else add to R
			for (int elem = 0; elem < load.apply_to.size() / 2; elem++)
				elements[load.apply_to[2 * elem] - 1]->set_load(load.type, load.apply_to[2 * elem + 1], load.data);
			
	}
	logger& log = logger::log();
	log.print("Create loads");
}

void Data::fillGlobalK() {
	int inf_count = 0; // ??
	int nodes_count = parser->mesh.nodes_count;
	int elems_count = parser->mesh.elems_count;

	K.resize(dim * (nodes_count + inf_count), dim * (nodes_count + inf_count));
	std::vector <Eigen::Triplet <double>> tripl_vec;
	for (int i = 0; i < elems_count; i++) {
		Eigen::MatrixXd loc_k = elements[i]->localK();
		for (int j = 0; j < elements[i]->nodes_count() * dim; j++)
			for (int k = 0; k < elements[i]->nodes_count() * dim; k++) {
				Eigen::Triplet <double> trpl(dim * (elements[i]->get_nodes(j / dim) - 1) + j % dim, dim * (elements[i]->get_nodes(k / dim) - 1) + k % dim, loc_k(j, k));
				tripl_vec.push_back(trpl);
			}
	}

	K.setFromTriplets(tripl_vec.begin(), tripl_vec.end());	
	zeroDiagonalCheck();
	logger& log = logger::log();
	log.print("Fill stiffness matrix");

} 

void Data::fillGlobalR() {
	int inf_count = 0; // ??
	int nodes_count = parser->mesh.nodes_count;
	int elems_count = parser->mesh.elems_count;

	R.resize(dim * (nodes_count + inf_count));
	for (int i = 0; i < elems_count; i++) {
		std::vector loc_r = elements[i]->localR();
		for (int j = 0; j < elements[i]->nodes_count() * dim; j++)
			R.insert(dim * (elements[i]->get_nodes(j / dim) - 1) + j % dim) = loc_r[j];
	}
	logger& log = logger::log();
	log.print("Fill right vector");
}

void Data::fillconstraints() {
	// c.first - comp, c.second - value 
	for (int node = 0; node < nodes.size(); node++)
		for (auto const& c: nodes[node].constraints)
			for (int i = 0; i < K.outerSize(); i++) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
					if (((it.row() == node * dim + c.first) ||
						(it.col() == node * dim + c.first)) && (it.row() != it.col())) {
						it.valueRef() = 0.0;
					}
					else if (((it.row() == node * dim + c.first) ||
						(it.col() == node * dim + c.first)) && (it.row() == it.col()))
						R.coeffRef(it.row()) = c.second * it.value(); 
			}
	logger& log = logger::log();
	log.print("Fill canctriants");
}

void Data::addToGlobalK(int first_index, int second_index, double value) {
	Eigen::Triplet <double> tripl (first_index, second_index, value);
	K.setFromTriplets(&tripl, &tripl + 1);
}

void Data::addToGlobalR(int index, double value) {
	R.insert(index) = value;
}

void Data::zeroDiagonalCheck() {
	for (int i = 0; i < nodes.size() * dim; i++)
		if (K.coeffRef(i, i) == 0) {
			std::cout << "Zero on diagonal: node " + std::to_string(static_cast<int>(i / dim + 1)) + " dof " + std::to_string(static_cast<int>(i % dim)) << std::endl;
			break;
		}
}

void Data::solve() {
	fillGlobalK();
	fillGlobalR();
	fillconstraints();

	zeroDiagonalCheck();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(K);

	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//solver.analyzePattern(K);
	//solver.factorize(K);
	
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	//solver.compute(K);

	if (solver.info() != Eigen::Success)
		std::cerr << "Error in K" << std::endl;

	U = solver.solve(R);

	logger& log = logger::log();
	log.print("Solve done");

}

void Data::writeMatrix(Eigen::SparseMatrix<double> A) {
	for (int i = 0; i < A.innerSize(); i++) {
		for (int j = 0; j < A.outerSize(); j++)
			std::cout << A.coeffRef(i, j) << " ";
		std::cout << std::endl;
	}
}

void Data::writeVector(Eigen::SparseVector<double> B) {
	for (int i = 0; i < B.size(); i++)
		std::cout << B.coeffRef(i) << std::endl;
}
