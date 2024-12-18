#include "Data.h"

Data::Data(std::shared_ptr <Parser> p) : parser(p) {
	this->dim = 2; // need read from fc
	logger& log = logger::log();
	log.print("Start parsing");

	create_nodes();
	create_elements();
	create_constants();
	create_constraints();
	create_load();

	log.print("End parsing");
}

void Data::set_output_param(std::vector<double> start, std::vector<double> end, int count) {
	this->line_end = end;
	this->line_start = start;
	this->points_count = count;
}

void Data::create_nodes() {
	for (int i = 0; i < parser->mesh.nodes_count; i++) {
		std::array <double, 3> coords;
		for (int j = 0; j < 3; j++)
			coords[j] = parser->mesh.nodes_coord[3 * i + j];
		this->nodes.push_back(Node(parser->mesh.node_id[i], coords));
	}
}

void Data::create_elements() {
	int node_tmp = 0;
	for (int i = 0; i < parser->mesh.elems_count; i++) {
		std::vector<int> elem_nodes;
		std::shared_ptr<Element> elem;
		switch (parser->mesh.elem_types[i]) {
		case TRI:
			elem_nodes.resize(3);
			for (int j = 0; j < 3; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 3;
			elem = std::make_shared<triElement>(triElement(parser->mesh.elem_id[i], TRI, elem_nodes));
			break;
		case QUAD:
			elem_nodes.resize(4);
			for (int j = 0; j < 4; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 4;
			elem = std::make_shared<quadElement>(quadElement(parser->mesh.elem_id[i], QUAD, elem_nodes));
			break;
		case TETRA:
			elem_nodes.resize(4);
			for (int j = 0; j < 4; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 4;
			elem = std::make_shared<tetraElement>(tetraElement(parser->mesh.elem_id[i], TETRA, elem_nodes));
			break;
		case HEX:
			elem_nodes.resize(8);
			for (int j = 0; j < 8; j++)
				elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			node_tmp += 8;
			elem = std::make_shared<hexElement>(hexElement(parser->mesh.elem_id[i], HEX, elem_nodes));
			break;
		default:
			throw runtime_error("Error: incorrect element " + to_string(i) +
				" type: " + to_string(parser->mesh.elem_types[i]));
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
}

void Data::create_constants() {
	// need block
	for (int i = 0; i < parser->mesh.elems_count; i++)
		elements[i]->set_constants(parser->material[0].constants[0], parser->material[0].constants[1]);
}

void Data::create_constraints() {
	for (int id = 0; id < parser->restraints.size(); id++) {
		auto& restraints = parser->restraints[id];
		for (int node = 0; node < restraints.size; node++)
			for (int i = 0; i < 6; i++)
				if (restraints.flag[i])
					nodes[restraints.apply_to[node] - 1].set_constraints(i, restraints.data[i]); // if exist --> throw
	}
}

void Data::create_load() { // type - pressure
	for (int id = 0; id < parser->load.size(); id++) {
		auto& load = parser->load[id];
		if (load.type == PRESSURE) 		// else add to F
			for (int elem = 0; elem < load.apply_to.size() / 2; elem++)
				elements[load.apply_to[2 * elem] - 1]->set_load(load.type, load.apply_to[2 * elem + 1], load.data);
			
	}
}

void Data::fillGlobalK() {
	logger& log = logger::log();
	log.print("Start filling stiffness matrix");
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
	log.print("End filling stiffness matrix");
} 

void Data::fillGlobalF() {
	logger& log = logger::log();
	log.print("Start filling right vector");

	int inf_count = 0; // ??
	int nodes_count = parser->mesh.nodes_count;
	int elems_count = parser->mesh.elems_count;

	F.resize(dim * (nodes_count + inf_count));
	for (int i = 0; i < elems_count; i++) {
		std::vector<double> loc_f = elements[i]->localF();
		for (int j = 0; j < elements[i]->nodes_count() * dim; j++)
			F.coeffRef(dim * (elements[i]->get_nodes(j / dim) - 1) + j % dim) += loc_f[j];
	}
	log.print("End filling right vector");
}

void Data::fillconstraints() {
	logger& log = logger::log();
	log.print("Start filling constraints");
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
						F.coeffRef(it.row()) = c.second * it.value(); 
			}
	log.print("End filling constraints");
}

void Data::fillGlobalC() {
	int inf_count = 0; // ??
	int nodes_count = parser->mesh.nodes_count;
	int elems_count = parser->mesh.elems_count;

	C.resize(nodes_count + inf_count, nodes_count + inf_count);
	std::vector <Eigen::Triplet <double>> tripl_vec;
	for (int i = 0; i < elems_count; i++) {
		Eigen::MatrixXd loc_c = elements[i]->localC();
		for (int j = 0; j < elements[i]->nodes_count(); j++)
			for (int k = 0; k < elements[i]->nodes_count(); k++) {
				Eigen::Triplet <double> trpl(elements[i]->get_nodes(j) - 1, elements[i]->get_nodes(k) - 1, loc_c(j, k));
				tripl_vec.push_back(trpl);
			}
	}

	C.setFromTriplets(tripl_vec.begin(), tripl_vec.end());
	zeroDiagonalCheck();
}

void Data::fillGlobalR(int type, int comp) {
	int inf_count = 0; // ??
	int nodes_count = parser->mesh.nodes_count;
	int elems_count = parser->mesh.elems_count;

	R.resize(nodes_count + inf_count);
	for (int i = 0; i < elems_count; i++) {
		std::vector<double> value;
		for (int j = 0; j < elements[i]->nodes_count(); j++)
			value.push_back(elements[i]->results[j][type](comp));
		std::vector<double> loc_r = elements[i]->localR(value);
		for (int j = 0; j < elements[i]->nodes_count(); j++)
			R.coeffRef(elements[i]->get_nodes(j) - 1) += loc_r[j];
	}
}

void Data::addToGlobalK(int first_index, int second_index, double value) {
	Eigen::Triplet <double> tripl (first_index, second_index, value);
	K.setFromTriplets(&tripl, &tripl + 1);
}

void Data::addToGlobalF(int index, double value) {
	F.coeffRef(index) = value;
}

void Data::displacementToElements() {
	std::vector <double> disp;
	for (int elem = 0; elem < elements.size(); elem++) {
		elements[elem]->results.resize(elements[elem]->nodes_count());
		for (int node = 0; node < elements[elem]->nodes_count(); node++) {
			elements[elem]->results[node].resize(COUNT);
			//elements[elem]->results[node][DISPLACEMENT].resize(dim);

			//elements[elem]->results[node][DISPLACEMENT][X] = U(dim * (elements[elem]->get_nodes(node) - 1));
			//elements[elem]->results[node][DISPLACEMENT][Y] = U(dim * (elements[elem]->get_nodes(node) - 1) + 1);
			//if (dim == 3)
			//	elements[elem]->results[node][DISPLACEMENT][Z] = U(dim * (elements[elem]->get_nodes(node) - 1) + 2);

			elements[elem]->displacements.resize(dim * elements[elem]->nodes_count());
			elements[elem]->displacements[dim * node] = U(dim * (elements[elem]->get_nodes(node) - 1));
			elements[elem]->displacements[dim * node + 1] = U(dim * (elements[elem]->get_nodes(node) - 1) + 1);
			if (dim == 3)
				elements[elem]->displacements[dim * node + 2] = U(dim * (elements[elem]->get_nodes(node) - 1) + 2);

		}
	}
}

void Data::displacementToNodes() {
	for (int i = 0; i < nodes.size(); i++)
		for (int j = 0; j < dim; j++) {
			nodes[i].set_res_size(DISPLACEMENT, dim);
			nodes[i].set_result(U.coeffRef(dim * i + j), DISPLACEMENT, j);
		}
}

void Data::calcStrain() {
	for (int elem = 0; elem < elements.size(); elem++) {
		//std::vector <double> ksi = { 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774, 0.5774 };
		//std::vector <double> eta = { 0.5774, 0.5774, -0.5774, -0.5774, 0.5774, 0.5774, -0.5774, -0.5774 };
		//std::vector <double> zeta = { 0.5774, 0.5774, 0.5774, 0.5774, -0.5774, -0.5774, -0.5774, -0.5774 };
		std::vector <double> ksi = { 1, -1, -1, 1, 1, -1, -1, 1 };
		std::vector <double> eta = { 1, 1, -1, -1, 1, 1, -1, -1 };
		std::vector <double> zeta = { 1, 1, 1, 1, -1, -1, -1, -1 };
		//int count;
		//if (elements[elem]->get_type() == TRI || TETRA)
		//	count == 1;
		//else
		//	count == elements[elem]->nodes_count();

		for (int node = 0; node < elements[elem]->nodes_count(); node++) {
			//elements[elem]->results[node][STRAIN].resize(output_fields(STRAIN, dim), 0);
			elements[elem]->results[node][STRAIN] = elements[elem]->B(ksi[node], eta[node], zeta[node]) * elements[elem]->displacements;
			//productMV(elements[elem]->B(ksi[node], eta[node], zeta[node]), elements[elem]->displacements, elements[elem]->results[node][STRAIN]);
		}
	}
	logger& log = logger::log();
	log.print("Calculate Strain");
}

void Data::calcStress() {
	for (int elem = 0; elem < elements.size(); elem++) {
		for (int node = 0; node < elements[elem]->nodes_count(); node++) {
			//elements[elem]->results[node][STRESS].resize(output_fields(STRESS, dim), 0);
			elements[elem]->results[node][STRESS] = elements[elem]->planeStrainD() * elements[elem]->results[node][STRAIN];
			if (dim == 2)
				elements[elem]->results[node][STRAIN][XY_2D] /= 2;

			//productMV(elements[elem]->planeStrainD(), elements[elem]->results[node][STRAIN], elements[elem]->results[node][STRESS]);
		}
	}
	logger& log = logger::log();
	log.print("Calculate Stress");
}

void Data::zeroDiagonalCheck() {
	for (int i = 0; i < nodes.size() * dim; i++)
		if (K.coeffRef(i, i) == 0) {
			std::cout << "Zero on diagonal: node " + std::to_string(static_cast<int>(i / dim + 1)) + " dof " + std::to_string(static_cast<int>(i % dim)) << std::endl;
			break;
		}
}

void Data::outputValues(int type, int comp) {
	std::vector<double> values;
	values.resize(points_count, 0);
	std::vector<std::vector<double>> points;
	output_points(line_start, line_end, points_count, points);
	interpolation(points, values, type, comp);

	ofstream file;
	std::string fieldname;
	field_name(fieldname, type, comp, dim);
	file.open(fieldname + ".txt");

	double line_len = sqrt(pow((line_end[0] - line_start[0]), 2) + pow((line_end[1] - line_start[1]), 2));
	for (int i = 0; i < points_count; i++) 
		//file << "(" << line_len / points_count * i << ", " << std::fixed << std::setprecision(12) << values[i] << ")" << std::endl;
		//file << points[i][0] << std::endl;
		file << std::fixed << std::setprecision(12) << values[i] << std::endl;
		
	file.close();

	if (type == STRESS)
		for (int i = 0; i < points_count; i++)
			out_stress[comp].push_back(values[i]);

	logger& log = logger::log();
	log.print("Interpolation " + fieldname);
}

void Data::interpolation(std::vector<std::vector<double>>& points, 
						std::vector<double>& values, int type, int comp) {
	for (int p = 0; p < points_count; p++)
		for (int i = 0; i < elements.size(); i++) {
			if (elements[i]->pointInElem(points[p])) {
				//values[p] = elements[i]->results[0][type][comp];

				std::vector<double> N = elements[i]->FF(points[p][0], points[p][1]);
				for (int node = 0; node < elements[i]->nodes_count(); node++)
					values[p] += N[node] * nodes[elements[i]->get_nodes(node) - 1].get_result(type, comp);
					//values[p] += N[node] * elements[i]->results[node][type][comp];

				break;
			}
		}
}

void Data::solve() {

	fillGlobalK();
	fillGlobalF();
	fillconstraints();

	zeroDiagonalCheck();

	logger& log = logger::log();
	log.print("Start solving");

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
	solver.compute(K);

	if (solver.info() != Eigen::Success)
		throw runtime_error("Error in K");

	U = solver.solve(F);
	log.print("Solve done");

	fillFields();
	smoothing();
}

void Data::fillFields() {
	displacementToElements();
	calcStrain();
	calcStress();
}

void Data::smoothing() {
	logger& log = logger::log();
	log.print("Start smoothing");

	fillGlobalC();
	displacementToNodes();

	out_stress.resize((dim == 2) ? 3 : 6);

	for (int type = 1; type < COUNT; type++)
		for (int comp = 0; comp < output_fields(type, dim); comp++) {
			fillGlobalR(type, comp);

			Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
			//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
			solver.compute(C);
			if (solver.info() != Eigen::Success)
				throw runtime_error("Error in C");

			Eigen::MatrixXd Result;
			Result = solver.solve(R);

			for (int i = 0; i < nodes.size(); i++) {
				nodes[i].set_res_size(type, dim);
				nodes[i].set_result(Result(i), type, comp);
			}

			std::string fieldname;
			field_name(fieldname, type, comp, dim);
			log.print("Solve agreed resultant " + fieldname + " done");

			outputValues(type, comp);
		}
	for (int i = 0; i < output_fields(DISPLACEMENT, dim); i++)
		outputValues(DISPLACEMENT, i);

	log.print("End smoothing");
}

