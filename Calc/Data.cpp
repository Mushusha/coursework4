#include "Data.h"

Data::Data(const Data& other)
	: dim(other.dim),
	filename(other.filename),
	analisys_type(other.analisys_type),
	damping(other.damping),
	max_time(other.max_time),
	max_iter(other.max_iter),
	iter_res_output(other.iter_res_output),
	omega(other.omega),
	Amp(other.Amp),
	elements(other.elements),
	nodes(other.nodes) {
}
Data& Data::operator=(const Data& other) {
	if (this != &other) {
		dim = other.dim;
		filename = other.filename;
		analisys_type = other.analisys_type;
		damping = other.damping;
		max_time = other.max_time;
		max_iter = other.max_iter;
		iter_res_output = other.iter_res_output;
		omega = other.omega;
		Amp = other.Amp;
		elements = other.elements;
		nodes = other.nodes;
	}
	return *this;
}
Data::Data(Data&& other) noexcept
	: dim(std::exchange(other.dim, 0)),
	filename(std::move(other.filename)),
	analisys_type(std::exchange(other.analisys_type, 0)),
	damping(std::exchange(other.damping, 0)),
	max_time(std::exchange(other.max_time, 0)),
	max_iter(std::exchange(other.max_iter, 0)),
	iter_res_output(std::exchange(other.iter_res_output, 0)),
	omega(std::exchange(other.omega, 0)),
	Amp(std::exchange(other.Amp, 0)),
	elements(std::move(other.elements)),
	nodes(std::move(other.nodes)) {
}
Data& Data::operator=(Data&& other) noexcept {
	if (this != &other) {
		dim = std::exchange(other.dim, 0);
		filename = std::move(other.filename);
		analisys_type = std::exchange(other.analisys_type, 0);
		damping = std::exchange(other.damping, 0);
		max_time = std::exchange(other.max_time, 0);
		max_iter = std::exchange(other.max_iter, 0);
		iter_res_output = std::exchange(other.iter_res_output, 0);
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

	filename = parser->get_filename();

	analisys_type = parser->settings.analisys_type;
	damping = parser->settings.d;
	max_time = parser->settings.max_time;
	max_iter = parser->settings.max_iter;
	iter_res_output = parser->settings.iter_res_output;
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
				elem_nodes.resize(count_nodes(TRI));
				for (int j = 0; j < elem_nodes.size(); j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				elem = std::make_shared<Tri>(Tri(parser->mesh.elem_id[i], TRI, elem_nodes));
				break;
			case QUAD:
				elem_nodes.resize(count_nodes(QUAD));
				for (int j = 0; j < elem_nodes.size(); j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				elem = std::make_shared<Quad>(Quad(parser->mesh.elem_id[i], QUAD, elem_nodes));
				break;
			case TETRA:
				elem_nodes.resize(count_nodes(TETRA));
				for (int j = 0; j < elem_nodes.size(); j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				elem = std::make_shared<Tetra>(Tetra(parser->mesh.elem_id[i], TETRA, elem_nodes));
				break;
			case HEX:
				elem_nodes.resize(count_nodes(HEX));
				for (int j = 0; j < elem_nodes.size(); j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				elem = std::make_shared<Hex>(Hex(parser->mesh.elem_id[i], HEX, elem_nodes));
				break;
			case WEDGE:
				elem_nodes.resize(count_nodes(WEDGE));
				for (int j = 0; j < elem_nodes.size(); j++)
					elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
				elem = std::make_shared<Wedge>(Wedge(parser->mesh.elem_id[i], WEDGE, elem_nodes));
				break;
			//case PYR:
			//	elem_nodes.resize(count_nodes(PYR));
			//	for (int j = 0; j < elem_nodes.size(); j++)
			//		elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
			//	elem = std::make_shared<Pyr>(Pyr(parser->mesh.elem_id[i], PYR, elem_nodes));
			//	break;
			default:
				throw runtime_error("Error: incorrect element " + to_string(i) +
									" type: " + to_string(parser->mesh.elem_types[i]));
				break;
		}
		node_tmp += elem_nodes.size();

		std::vector<double> x, y, z;
		for (int j = 0; j < elem_nodes.size(); j++) {
			x.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 0]);
			y.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 1]);
			z.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 2]); // may be Nodes 
		}
		elem->set_coords(x, y, z);
		this->elements.push_back(elem);
	}
	create_infelements(parser);
}

void Data::create_infelements(std::shared_ptr <const Parser> parser) {
	num_inf_elems = 0;

	if (parser->infinite.size() == 0)
		return;

	const auto inf = parser->infinite[0];

	std::map<int, int> map_node_inf;
	std::map<int, std::array<int, 3>> node_coords; // [i] - num node, [j] - x | y | z

	for (int i = 0; i < inf.size; i++) {
		const int elem = inf.apply_to[2 * i];
		const int side = inf.apply_to[2 * i + 1];

		int num_nodes;
		std::vector<int> side_nodes = this->elements[elem - 1]->edge_to_node(side);
		std::vector<int> new_nodes;

		for (int& node : side_nodes) {
			int glob_node = this->elements[elem - 1]->get_node(node);

			if (map_node_inf.find(glob_node) == map_node_inf.end()) {
				map_node_inf[glob_node] = nodes.size() + 1;
				std::array <double, 3> coords;
				for (int i = 0; i < 3; i++)
					coords[i] = 2 * this->nodes[glob_node - 1]->getCoord(i) + this->nodes[inf.point - 1]->getCoord(i);

				this->nodes.push_back(std::make_shared<Node>(map_node_inf.at(glob_node), coords));
			}

			new_nodes.push_back(glob_node);
			new_nodes.push_back(map_node_inf.find(glob_node)->second);
		}

		// ONLY INFQUAD
		std::vector<int> inf_elem_nodes = { new_nodes[0], new_nodes[1], new_nodes[3], new_nodes[2] };

		std::vector<double> x1, y1, z1;
		for (auto& node : inf_elem_nodes) {
			x1.push_back(this->nodes[node - 1]->getX());
			y1.push_back(this->nodes[node - 1]->getY());
			z1.push_back(this->nodes[node - 1]->getZ());
		}

		std::shared_ptr<Element> new_elem = std::make_shared<infQuad>(infQuad(this->elements.size(), INFQUAD, inf_elem_nodes));
		new_elem->set_coords(x1, y1, z1);
		this->elements.push_back(new_elem);

		if (parser->settings.analisys_type == "dynamic") {
			num_inf_elems++; // complex
			dynamic_cast<infQuad*>(new_elem.get())->is_dyn = true;
			dynamic_cast<infQuad*>(new_elem.get())->omega = omega;
		}

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
			for (int i = 0; i < 3; i++)
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

			double factor = Young * (1 - Poisson) / ((1 + Poisson) * (1 - 2 * Poisson));

			this->elements[i]->D = Eigen::MatrixXd::Zero(6, 6);

			this->elements[i]->D(0, 0) = 1.0;
			this->elements[i]->D(0, 1) = Poisson / (1 - Poisson);
			this->elements[i]->D(0, 2) = Poisson / (1 - Poisson);

			this->elements[i]->D(1, 1) = 1.0;
			this->elements[i]->D(1, 2) = Poisson / (1 - Poisson);

			this->elements[i]->D(2, 2) = 1.0;

			this->elements[i]->D(3, 3) = (1 - 2 * Poisson) / (2 * (1 - Poisson));
			this->elements[i]->D(4, 4) = (1 - 2 * Poisson) / (2 * (1 - Poisson));
			this->elements[i]->D(5, 5) = (1 - 2 * Poisson) / (2 * (1 - Poisson));

			this->elements[i]->D *= factor;

			for (int j = 0; j < 6; j++)
				for (int k = j + 1; k < 6; k++)
					this->elements[i]->D(k, j) = this->elements[i]->D(j, k);
		}
}