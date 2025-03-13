#include "Data.h"

Data::Data(const Data& other)
	: dim(other.dim),
	line_start(other.line_start),
	line_end(other.line_end),
	points_count(other.points_count),
	out_stress(other.out_stress),
	elements(other.elements),
	nodes(other.nodes),
	parser(other.parser) {
}
Data& Data::operator=(const Data& other) {
	if (this != &other) {
		dim = other.dim;
		line_start = other.line_start;
		line_end = other.line_end;
		points_count = other.points_count;
		out_stress = other.out_stress;
		elements = other.elements;
		nodes = other.nodes;
		parser = other.parser;
	}
	return *this;
}
Data::Data(Data&& other) noexcept
	: dim(std::exchange(other.dim, 0)),
	line_start(std::move(other.line_start)),
	line_end(std::move(other.line_end)),
	points_count(std::exchange(other.points_count, 0)),
	out_stress(std::move(other.out_stress)),
	elements(std::move(other.elements)),
	nodes(std::move(other.nodes)),
	parser(std::move(other.parser)) {
}
Data& Data::operator=(Data&& other) noexcept {
	if (this != &other) {
		dim = std::exchange(other.dim, 0);
		line_start = std::move(other.line_start);
		line_end = std::move(other.line_end);
		points_count = std::exchange(other.points_count, 0);
		out_stress = std::move(other.out_stress);
		elements = std::move(other.elements);
		nodes = std::move(other.nodes);
		parser = std::move(other.parser);
	}
	return *this;
}

Data::Data(std::shared_ptr <Parser> p) : parser(p) {
	logger& log = logger::log();
	log.print("Start parsing");

	if (parser->settings.dimensions == "2D")
		this->dim = 2;
	else
		this->dim = 3;

	analisys_type = parser->settings.analisys_type;

	create_nodes();
	create_elements();
	create_constants();
	create_constraints();
	create_load();
	create_D();

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
	create_infelements();
	parser->mesh.elems_count = this->elements.size();
	parser->mesh.nodes_count = this->nodes.size();
}

void Data::create_infelements() {
	// if 2D
	int inf_node = 0;
	for (int i = 0; i < parser->sidesets.size(); i++) {
		for (int j = 0; j < parser->sidesets[i].size; j++) {
			auto inf = parser->sidesets[i].apply_to;

			std::pair<double, double> C1(this->nodes[parser->nodesets[0].apply_to[inf_node]].getX(), this->nodes[parser->nodesets[0].apply_to[inf_node]].getY()); // nodesets[0] ??? 
			std::pair<double, double> C2(this->nodes[parser->nodesets[0].apply_to[inf_node + 1]].getX(), this->nodes[parser->nodesets[0].apply_to[inf_node + 1]].getY());

			std::pair<int, int> edge;
			if (inf[2 * j + 1] != 3)
				edge = std::make_pair(inf[2 * j + 1], inf[2 * j + 1] + 1);
			else if (inf[2 * j + 1] == 3)
				edge = std::make_pair(inf[2 * j + 1], 0);

			edge.first = this->elements[inf[2 * j]]->get_nodes(edge.first);
			edge.second = this->elements[inf[2 * j]]->get_nodes(edge.second);

			int n = this->nodes.size();
			std::vector<int> elem_nodes = { static_cast<int>(parser->nodesets[0].apply_to[inf_node]), n + 1, n + 2, static_cast<int>(parser->nodesets[0].apply_to[inf_node] + 1) };
			std::vector<double> x, y, z = { 0, 0, 0, 0 };

			x.push_back(this->nodes[edge.first - 1].getX() * 2 - C1.first);
			x.push_back(C1.first);
			x.push_back(C2.first);
			x.push_back(this->nodes[edge.second - 1].getX() * 2 - C2.first);

			y.push_back(this->nodes[edge.first - 1].getY() * 2 - C1.second);
			y.push_back(C1.second);
			y.push_back(C2.second);
			y.push_back(this->nodes[edge.second - 1].getY() * 2 - C2.second);

			std::shared_ptr<Element> elem = std::make_shared<infQuadElement>(infQuadElement(this->elements.size(), INFQUAD, elem_nodes));
			elem->set_coords(x, y, z);
			this->elements.push_back(elem);

			inf_node++;
			for (int k = 2; k < 4; k++) {
				std::array <double, 3> coords = { x[k], y[k], z[k] };
				this->nodes.push_back(Node(elem_nodes[k], coords));
			}

			if (parser->sidesets[i].load != -1) { // type == PRESSURE
				auto& load = parser->load[parser->sidesets[i].load];
				elements[this->elements.size() - 1]->set_load(load.type, 3, load.data);
			}
		}
	}
}

void Data::create_constants() {
	// need block
	for (int i = 0; i < parser->mesh.elems_count; i++)
		elements[i]->set_constants(parser->material[0].constants[0], parser->material[0].constants[1], parser->material[0].constants[2]);
}

void Data::create_constraints() {
	for (int id = 0; id < parser->restraints.size(); id++) {
		auto& restraints = parser->restraints[id];
		for (int node = 0; node < restraints.size; node++)
			for (int i = 0; i < 6; i++)
				if (restraints.flag[i])
					nodes[restraints.apply_to[node] - 1].set_constraints(i, restraints.data[i]);
	}
}

void Data::create_load() { // type - pressure
	for (int id = 0; id < parser->load.size(); id++)
		if (!parser->load[id].inf) {
			auto& load = parser->load[id];
			if (load.type == PRESSURE) 		// else add to F
				for (int elem = 0; elem < load.apply_to.size() / 2; elem++)
					elements[load.apply_to[2 * elem] - 1]->set_load(load.type, load.apply_to[2 * elem + 1], load.data);
			else if (load.type == NODEFORCE)
				for (int i = 0; i < 6; i++)
					if (load.data[i] == 0.0)
						continue;
					else
						for (auto node : load.apply_to)
							nodes[node].load[i] = load.data[i];
		}
}

void Data::create_D() {
	if (parser->settings.dimensions == "2D")
		if (parser->settings.plane_state == "p-stress")
			for (int i = 0; i < parser->mesh.elems_count; i++) {
				double Poisson = this->elements[i]->get_nu();
				double Young = this->elements[i]->get_E();
				this->elements[i]->D = Eigen::MatrixXd::Zero(3, 3);
				this->elements[i]->D(0, 0) = 1;
				this->elements[i]->D(0, 1) = Poisson;
				this->elements[i]->D(1, 0) = Poisson;
				this->elements[i]->D(1, 1) = 1;
				this->elements[i]->D(2, 2) = (1 - Poisson) / 2;

				this->elements[i]->D *= Young / (1 - pow(Poisson, 2));
			}
		else
			for (int i = 0; i < parser->mesh.elems_count; i++) {
				double Poisson = this->elements[i]->get_nu();
				double Young = this->elements[i]->get_E();

				this->elements[i]->D = Eigen::MatrixXd::Zero(3, 3);
				this->elements[i]->D(0, 0) = 1;
				this->elements[i]->D(0, 1) = Poisson / (1 - Poisson);
				this->elements[i]->D(1, 0) = Poisson / (1 - Poisson);
				this->elements[i]->D(1, 1) = 1;
				this->elements[i]->D(2, 2) = (1 - 2 * Poisson) / (2 * (1 - Poisson));

				this->elements[i]->D *= Young * (1 - Poisson) / ((1 + Poisson) * (1 - 2 * Poisson));
			}
	else
		for (int i = 0; i < parser->mesh.elems_count; i++) {
			double Poisson = this->elements[i]->get_nu();
			double Young = this->elements[i]->get_E();

			this->elements[i]->D = Eigen::MatrixXd::Zero(6, 6);
			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 3; j++) {
					if (i == j && i < 3)
						this->elements[i]->D(i, j) = 1;
					if (i != j && i < 3)
						this->elements[i]->D(i, j) = Poisson / (1 - Poisson);
					if (i >= 3)
						this->elements[i]->D(i, i) = (1 - 2 * Poisson) / (2 * (1 - Poisson));
				}

			this->elements[i]->D *= Young * (1 - Poisson) / ((1 + Poisson) * (1 - 2 * Poisson));
		}
}

//void Data::outputValues(int type, int comp) {
//	std::vector<double> values;
//	values.resize(points_count, 0);
//	std::vector<std::vector<double>> points;
//	output_points(line_start, line_end, points_count, points);
//	interpolation(points, values, type, comp);
//
//	ofstream file;
//	std::string fieldname;
//	field_name(fieldname, type, comp, dim);
//	file.open(fieldname + ".txt");
//
//	double line_len = sqrt(pow((line_end[0] - line_start[0]), 2) + pow((line_end[1] - line_start[1]), 2));
//	for (int i = 0; i < points_count; i++) 
//		//file << "(" << line_len / points_count * i << ", " << std::fixed << std::setprecision(12) << values[i] << ")" << std::endl;
//		file << "(" << points[i][0] << ", " << std::fixed << std::setprecision(12) << values[i] << ")" << std::endl;
//		//file << points[i][0] << std::endl;
//		//file << std::fixed << std::setprecision(12) << values[i] << std::endl;
//		
//	file.close();
//
//	if (type == STRESS)
//		for (int i = 0; i < points_count; i++)
//			out_stress[comp].push_back(values[i]);
//
//	logger& log = logger::log();
//	log.print("Interpolation " + fieldname);
//}

//void Data::printMeshStress() {
//	std::ofstream file_coord_x;
//	file_coord_x.open("coord_x.txt");
//	std::ofstream file_coord_y;
//	file_coord_y.open("coord_y.txt");
//	//std::ofstream file_coord_z;
//	//file_coord_z.open("coord_z.txt");
//
//	std::ofstream file_stress_xx;
//	file_stress_xx.open("all_stress_xx.txt");
//	std::ofstream file_stress_yy;
//	file_stress_yy.open("all_stress_yy.txt");
//	std::ofstream file_stress_xy;
//	file_stress_xy.open("all_stress_xy.txt");
//	//std::ofstream file_stress_xz;
//	//file_stress_xz.open("all_stress_xz.txt");
//	//std::ofstream file_stress_yz;
//	//file_stress_yz.open("all_stress_yz.txt");
//	//std::ofstream file_stress_zz;
//	//file_stress_zz.open("all_stress_zz.txt");
//
//	for (int i = 0; i < nodes.size(); i++)
//		nodes[i].printStress();
//
//	file_coord_x.close();
//	file_coord_y.close();
//	//file_coord_z.close();
//
//	file_stress_xx.close();
//	file_stress_yy.close();
//	file_stress_xy.close();
//	//file_stress_xz.close();
//	//file_stress_yz.close();
//	//file_stress_zz.close();
//
//	std::ofstream file_disp_x;
//	file_disp_x.open("all_disp_x.txt");
//	std::ofstream file_disp_y;
//	file_disp_y.open("all_disp_y.txt");
//
//	for (int i = 0; i < nodes.size(); i++) {
//		file_disp_x << U(dim * i) << std::endl;
//		file_disp_y << U(dim * i + 1) << std::endl;
//	}
//
//	file_disp_x.close();
//	file_disp_y.close();
//}

//void Data::interpolation(std::vector<std::vector<double>>& points, 
//						std::vector<double>& values, int type, int comp) {
//	for (int p = 0; p < points_count; p++)
//		for (int i = 0; i < elements.size(); i++) {
//			if (elements[i]->pointInElem(points[p])) {
//				//values[p] = elements[i]->results[0][type][comp];
//				std::vector<double> Coord = elements[i]->coordFF(points[p][0], points[p][1]); // dim == 3
//				std::vector<double> N = elements[i]->FF(Coord[0], Coord[1]);
//				for (int node = 0; node < elements[i]->nodes_count(); node++)
//					values[p] += N[node] * nodes[elements[i]->get_nodes(node) - 1].get_result(type, comp);
//
//					//if (type != DISPLACEMENT)
//					//	values[p] += N[node] * elements[i]->results[node][type](comp);
//					//else
//					//	values[p] += N[node] * elements[i]->displacements(node * dim + comp);
//				break;
//			}
//		}
//}

//void Data::smoothing() {
//	logger& log = logger::log();
//	log.print("Start smoothing");
//
//	fillGlobalC();
//	displacementToNodes();
//
//	out_stress.resize((dim == 2) ? 3 : 6);
//
//	for (int i = 0; i < output_fields(DISPLACEMENT, dim); i++)
//		outputValues(DISPLACEMENT, i);
//
//	for (int type = 1; type < COUNT; type++)
//		for (int comp = 0; comp < output_fields(type, dim); comp++) {
//			fillGlobalR(type, comp);
//
//			Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
//			//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
//			solver.compute(C);
//			if (solver.info() != Eigen::Success)
//				throw runtime_error("Error in C");
//
//			Eigen::MatrixXd Result;
//			Result = solver.solve(R);
//
//			for (int i = 0; i < nodes.size(); i++) {
//				nodes[i].set_res_size(type, dim);
//				nodes[i].set_result(Result(i), type, comp);
//			}
//
//			std::string fieldname;
//			field_name(fieldname, type, comp, dim);
//			log.print("Solve agreed resultant " + fieldname + " done");
//
//			outputValues(type, comp);
//		}
//
//	printMeshStress();
//
//	log.print("End smoothing");
//}

