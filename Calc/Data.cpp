#include "Data.h"

Data::Data(const Data& other)
	: dim(other.dim),
	out_stress(other.out_stress),
	analisys_type(other.analisys_type),
	damping(other.damping),
	max_time(other.max_time),
	max_iter(other.max_iter),
	omega(other.omega),
	Amp(other.Amp),
	elements(other.elements),
	nodes(other.nodes) {
}
Data& Data::operator=(const Data& other) {
	if (this != &other) {
		dim = other.dim;
		out_stress = other.out_stress;
		analisys_type = other.analisys_type;
		damping = other.damping;
		max_time = other.max_time;
		max_iter = other.max_iter;
		omega = other.omega;
		Amp = other.Amp;
		elements = other.elements;
		nodes = other.nodes;
	}
	return *this;
}
Data::Data(Data&& other) noexcept
	: dim(std::exchange(other.dim, 0)),
	out_stress(std::move(other.out_stress)),
	analisys_type(std::exchange(other.analisys_type, 0)),
	damping(std::exchange(other.damping, 0)),
	max_time(std::exchange(other.max_time, 0)),
	max_iter(std::exchange(other.max_iter, 0)),
	omega(std::exchange(other.omega, 0)),
	Amp(std::exchange(other.Amp, 0)),
	elements(std::move(other.elements)),
	nodes(std::move(other.nodes)) {
}
Data& Data::operator=(Data&& other) noexcept {
	if (this != &other) {
		dim = std::exchange(other.dim, 0);
		out_stress = std::move(other.out_stress);
		analisys_type = std::exchange(other.analisys_type, 0);
		damping = std::exchange(other.damping, 0);
		max_time = std::exchange(other.max_time, 0);
		max_iter = std::exchange(other.max_iter, 0);
		omega = std::exchange(other.omega, 0);
		Amp = std::exchange(other.Amp, 0);
		elements = std::move(other.elements);
		nodes = std::move(other.nodes);
	}
	return *this;
}

Data::Data(std::shared_ptr <const Parser> parser) {
	logger& log = logger::log();
	log.print("Start parsing");

	if (parser->settings.dimensions == "2D")
		this->dim = 2;
	else
		this->dim = 3;

	analisys_type = parser->settings.analisys_type;
	damping = parser->settings.d;
	max_time = parser->settings.max_time;
	max_iter = parser->settings.max_iter;
	omega = 10; // parser->
	Amp = 1e+08;

	create_nodes(parser);
	create_elements(parser);
	create_constants(parser);
	create_constraints(parser);
	create_load(parser);
	create_D(parser);

	log.print("End parsing");
}

void Data::create_nodes(std::shared_ptr <const Parser> parser) {
	for (int i = 0; i < parser->mesh.nodes_count; i++) {
		std::array <double, 3> coords;
		for (int j = 0; j < 3; j++)
			coords[j] = parser->mesh.nodes_coord[3 * i + j];
		this->nodes.push_back(std::make_shared<Node>(parser->mesh.node_id[i], coords));
	}
}

void Data::create_elements(std::shared_ptr <const Parser> parser) {
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
				elem = std::make_shared<Tri>(Tri(parser->mesh.elem_id[i], TRI, elem_nodes));
				break;
			case QUAD:
				elem_nodes.resize(4);
				for (int j = 0; j < 4; j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				node_tmp += 4;
				elem = std::make_shared<Quad>(Quad(parser->mesh.elem_id[i], QUAD, elem_nodes));
				break;
			case TETRA:
				elem_nodes.resize(4);
				for (int j = 0; j < 4; j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				node_tmp += 4;
				elem = std::make_shared<Tetra>(Tetra(parser->mesh.elem_id[i], TETRA, elem_nodes));
				break;
			case HEX:
				elem_nodes.resize(8);
				for (int j = 0; j < 8; j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				node_tmp += 8;
				elem = std::make_shared<Hex>(Hex(parser->mesh.elem_id[i], HEX, elem_nodes));
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
	create_infelements(parser);
}

void Data::create_infelements(std::shared_ptr <const Parser> parser) {
	// if 2D

	//int inf_node = 0;
	//for (int i = 0; i < parser->sidesets.size(); i++) {
	//	for (int j = 0; j < parser->sidesets[i].size; j++) {
	//		auto inf = parser->sidesets[i].apply_to;
	//
	//		//std::pair<double, double> C1(this->nodes[parser->nodesets[0].apply_to[inf_node]]->getX(), this->nodes[parser->nodesets[0].apply_to[inf_node]]->getY()); // nodesets[0] ??? 
	//		//std::pair<double, double> C2(this->nodes[parser->nodesets[0].apply_to[inf_node + 1]]->getX(), this->nodes[parser->nodesets[0].apply_to[inf_node + 1]]->getY());
	//
	//		std::pair<int, int> edge;
	//		if (inf[2 * j + 1] != this->elements[inf[2 * j]]->nodes_count() - 1)
	//			edge = std::make_pair(inf[2 * j + 1], inf[2 * j + 1] + 1);
	//		else
	//			edge = std::make_pair(inf[2 * j + 1], 0);
	//
	//		edge.first = this->elements[inf[2 * j]]->get_nodes(edge.first);
	//		edge.second = this->elements[inf[2 * j]]->get_nodes(edge.second);
	//
	//		int n = this->nodes.size();
	//		std::vector<int> elem_nodes = { static_cast<int>(parser->nodesets[0].apply_to[inf_node]), n + 1, n + 2, static_cast<int>(parser->nodesets[0].apply_to[inf_node] + 1) };
	//
	//		std::vector<double> x, y, z = { 0, 0, 0, 0 };
	//
	//		x.push_back(this->nodes[edge.first - 1]->getX() * 2 - C1.first);
	//		x.push_back(C1.first);
	//		x.push_back(C2.first);
	//		x.push_back(this->nodes[edge.second - 1]->getX() * 2 - C2.first);
	//
	//		y.push_back(this->nodes[edge.first - 1]->getY() * 2 - C1.second);
	//		y.push_back(C1.second);
	//		y.push_back(C2.second);
	//		y.push_back(this->nodes[edge.second - 1]->getY() * 2 - C2.second);
	//
	//		std::shared_ptr<Element> elem = std::make_shared<infQuad>(infQuad(this->elements.size(), INFQUAD, elem_nodes));
	//		elem->set_coords(x, y, z);
	//		this->elements.push_back(elem);
	//
	//		if (parser->settings.analisys_type == "dynamic") {
	//			dynamic_cast<infQuad*>(elem.get())->is_dyn = true;
	//			dynamic_cast<infQuad*>(elem.get())->omega = 10;
	//		}
	//		inf_node++;
	//		for (int k = 2; k < 4; k++) {
	//			std::array <double, 3> coords = { x[k], y[k], z[k] };
	//			this->nodes.push_back(make_shared<Node>(elem_nodes[k], coords));
	//		}
	//
	//		if (parser->sidesets[i].load != -1) { // type == PRESSURE
	//			auto& load = parser->load[parser->sidesets[i].load];
	//			elements[this->elements.size() - 1]->set_load(load.type, 3, load.data);
	//		}
	//	}
	//}

	num_inf_elems = 0;
	if (parser->nodesets.size() == 0)
		return;

	int n = this->nodes.size() + 1;
	// nodeset 0 - C, C1
	// nodeset 1 - Q, Q1
	for (int i = 0; i < parser->nodesets[0].apply_to.size() - 1; i++) {
		std::vector<int> elem_nodes = { n, static_cast<int>(parser->nodesets[1].apply_to[i]), static_cast<int>(parser->nodesets[1].apply_to[i + 1]), n + 1 };
		std::vector<double> x1, y1, z1 = { 0, 0, 0, 0 };

		x1.push_back(2 * this->nodes[static_cast<int>(parser->nodesets[0].apply_to[i] - 1)]->getX() - this->nodes[elem_nodes[1]]->getX());
		x1.push_back(this->nodes[elem_nodes[1] - 1]->getX());
		x1.push_back(this->nodes[elem_nodes[2] - 1]->getX());
		x1.push_back(2 * this->nodes[static_cast<int>(parser->nodesets[0].apply_to[i + 1] - 1)]->getX() - this->nodes[elem_nodes[2]]->getX());

		y1.push_back(2 * this->nodes[static_cast<int>(parser->nodesets[0].apply_to[i] - 1)]->getY() - this->nodes[elem_nodes[1]]->getY());
		y1.push_back(this->nodes[elem_nodes[1] - 1]->getY());
		y1.push_back(this->nodes[elem_nodes[2] - 1]->getY());
		y1.push_back(2 * this->nodes[static_cast<int>(parser->nodesets[0].apply_to[i + 1] - 1)]->getY() - this->nodes[elem_nodes[2]]->getY());

		std::shared_ptr<Element> elem = std::make_shared<infQuad>(infQuad(this->elements.size(), INFQUAD, elem_nodes));
		elem->set_coords(x1, y1, z1);
		this->elements.push_back(elem);

		std::array <double, 3> coords0 = { x1[0], y1[0], z1[0] };
		this->nodes.push_back(std::make_shared<Node>(elem_nodes[0], coords0));
		if (i == parser->nodesets[0].apply_to.size() - 2) {
			std::array <double, 3> coords3 = { x1[3], y1[3], z1[3] };
			this->nodes.push_back(std::make_shared<Node>(elem_nodes[3], coords3));
		}
		if (parser->settings.analisys_type == "dynamic") {
			num_inf_elems++; // complex
			dynamic_cast<infQuad*>(elem.get())->is_dyn = true;
			dynamic_cast<infQuad*>(elem.get())->omega = omega;
		}
		n++;
		num_inf_elems++;
	}
}

void Data::create_constants(std::shared_ptr <const Parser> parser) {
	// need block
	for (int i = 0; i < this->elements.size(); i++)
		elements[i]->set_constants(parser->material[0].constants[0], parser->material[0].constants[1], parser->material[0].constants[2]);
}

void Data::create_constraints(std::shared_ptr <const Parser> parser) {
	for (int id = 0; id < parser->restraints.size(); id++) {
		auto& restraints = parser->restraints[id];
		for (int node = 0; node < restraints.size; node++)
			for (int i = 0; i < 6; i++)
				if (restraints.flag[i])
					nodes[restraints.apply_to[node] - 1]->set_constraints(i, restraints.data[i]);
	}
}

void Data::create_load(std::shared_ptr <const Parser> parser) { // type - pressure
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
						for (int node = 0; node < load.apply_to.size(); node += 2)
							nodes[load.apply_to[node] - 1]->load[i] = load.data[i];
		}
}

void Data::create_D(std::shared_ptr <const Parser> parser) {
	if (parser->settings.dimensions == "2D")
		if (parser->settings.plane_state == "p-stress")
			for (int i = 0; i < this->elements.size(); i++) {
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
			for (int i = 0; i < this->elements.size(); i++) {
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
		for (int i = 0; i < this->elements.size(); i++) {
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