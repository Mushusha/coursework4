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
	create_constraints(parser);
	create_elements(parser);
	create_constants(parser);
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

void Data::create_elements(std::shared_ptr<const Parser> parser) {
    int node_tmp = 0;
    for (int i = 0; i < parser->mesh.elems_count; i++) {
        std::vector<int> elem_nodes;
        std::shared_ptr<Element> elem;

        const int order = parser->mesh.elem_orders[i];
        const ElemType type = static_cast<ElemType>(parser->mesh.elem_types[i]);

        elem_nodes.resize(count_nodes(type));
        for (int j = 0; j < elem_nodes.size(); j++)
            elem_nodes[j] = parser->mesh.elem_nodes[node_tmp + j];
        node_tmp += elem_nodes.size();

        if (type == QUADSEM && order > 2) {
            if (order == 3)
                elem = std::make_shared<SpectralQuad<4>>(SpectralQuad<4>(parser->mesh.elem_id[i], QUADSEM, elem_nodes));
            else if (order == 4)
                elem = std::make_shared<SpectralQuad<5>>(SpectralQuad<5>(parser->mesh.elem_id[i], QUADSEM, elem_nodes));
            else if (order == 5)
                elem = std::make_shared<SpectralQuad<6>>(SpectralQuad<6>(parser->mesh.elem_id[i], QUADSEM, elem_nodes));
            else
                throw std::runtime_error("Unsupported spectral element order: " + std::to_string(order));
            
			elem->set_order(order);

			build_spectral_quad_element(elem);
        }
		if (type == HEXSEM && order > 2) {
			if (order == 3)
				elem = std::make_shared<SpectralHex<4>>(SpectralHex<4>(parser->mesh.elem_id[i], HEXSEM, elem_nodes));
			else if (order == 4)
				elem = std::make_shared<SpectralHex<5>>(SpectralHex<5>(parser->mesh.elem_id[i], HEXSEM, elem_nodes));
			else if (order == 5)
				elem = std::make_shared<SpectralHex<6>>(SpectralHex<6>(parser->mesh.elem_id[i], HEXSEM, elem_nodes));
			else
				throw std::runtime_error("Unsupported spectral element order: " + std::to_string(order));

			elem->set_order(order);

			build_spectral_hex_element(elem);
		}
		else {
			switch (type) {
			case TRI:
				elem = std::make_shared<Tri>(Tri(parser->mesh.elem_id[i], TRI, elem_nodes));
				break;
			case QUAD:
				elem = std::make_shared<Quad>(Quad(parser->mesh.elem_id[i], QUAD, elem_nodes));
				break;
			case TETRA:
				elem = std::make_shared<Tetra>(Tetra(parser->mesh.elem_id[i], TETRA, elem_nodes));
				break;
			case HEX:
				elem = std::make_shared<Hex>(Hex(parser->mesh.elem_id[i], HEX, elem_nodes));
				break;
			case WEDGE:
				elem = std::make_shared<Wedge>(Wedge(parser->mesh.elem_id[i], WEDGE, elem_nodes));
				break;
			default:
				throw std::runtime_error("Error: incorrect element " + std::to_string(i) +
					" type: " + std::to_string(type));
			}

			std::vector<double> x, y, z;
			for (int j = 0; j < elem_nodes.size(); j++) {
				x.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 0]);
				y.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 1]);
				z.push_back(parser->mesh.nodes_coord[3 * (elem_nodes[j] - 1) + 2]);
			}

			elem->set_order(1);
			elem->set_coords(x, y, z);
		}
        this->elements.push_back(elem);
    }

    create_infelements(parser);
	renumber_nodes();
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
			for (int i = 0; i < dim; i++)
				if (restraints.flag[i])
					nodes[restraints.apply_to[node] - 1]->set_constraints(i, restraints.data[i]);
	}
}

void Data::create_load(std::shared_ptr <const Parser> parser) { // type - pressure
	for (int id = 0; id < parser->load.size(); id++)
		if (!parser->load[id].inf) {
			auto& load = parser->load[id];
			if (load.type == PRESSURE || load.type == PRESSURESEM) 		// else add to F
				for (int elem = 0; elem < load.apply_to.size() / 2; elem++)
					elements[load.apply_to[2 * elem] - 1]->set_load(load.type, load.apply_to[2 * elem + 1], load.data);
			else if (load.type == NODEFORCE) // add new spectral points
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

void Data::build_spectral_quad_element(std::shared_ptr<Element> elem) {
	int order = elem->get_order();
	int nodes_per_side = order + 1;

	std::vector<int> original_nodes = elem->get_node();
	if (original_nodes.size() != 8) {
		throw std::runtime_error("Expected QUAD8 element with 8 nodes for spectral element.");
	}

	std::map<int, std::map<int, double>> original_constraints;
	for (int node_id : original_nodes) {
		auto it = std::find_if(nodes.begin(), nodes.end(),
			[node_id](const std::shared_ptr<Node>& n) { return n->getID() == node_id; });
		if (it != nodes.end()) {
			original_constraints[node_id] = (*it)->constraints;
		}
	}

	std::vector<double> orig_x, orig_y, orig_z;
	for (int node_id : original_nodes) {
		auto it = std::find_if(nodes.begin(), nodes.end(),
			[node_id](const std::shared_ptr<Node>& n) { return n->getID() == node_id; });
		if (it == nodes.end()) {
			throw std::runtime_error("Node " + std::to_string(node_id) + " not found");
		}
		orig_x.push_back((*it)->getX());
		orig_y.push_back((*it)->getY());
		orig_z.push_back((*it)->getZ());
	}

	std::vector<double> gll_points, gll_weights;

	compute_gll_nodes_weights(order, gll_points, gll_weights);

	auto quad8_shape_functions = [](double xi, double eta, std::vector<double>& N) {
		N[0] = 0.25 * (1.0 - xi) * (1.0 - eta) * (-1.0 - xi - eta);
		N[1] = 0.25 * (1.0 + xi) * (1.0 - eta) * (-1.0 + xi - eta);
		N[2] = 0.25 * (1.0 + xi) * (1.0 + eta) * (-1.0 + xi + eta);
		N[3] = 0.25 * (1.0 - xi) * (1.0 + eta) * (-1.0 - xi + eta);

		N[4] = 0.5 * (1.0 - xi * xi) * (1.0 - eta);
		N[5] = 0.5 * (1.0 + xi) * (1.0 - eta * eta);
		N[6] = 0.5 * (1.0 - xi * xi) * (1.0 + eta);
		N[7] = 0.5 * (1.0 - xi) * (1.0 - eta * eta);
		};

	std::vector<double> N(8);
	std::vector<int> new_node_ids;
	std::vector<double> new_x, new_y, new_z;

	int next_node_id = nodes.size() + 1;
	if (!nodes.empty()) {
		int max_id = 0;
		for (const auto& node : nodes) {
			max_id = std::max(max_id, node->getID());
		}
		next_node_id = max_id + 1;
	}

	std::vector<bool> edge_fully_constrained(4, true);
	std::vector<std::vector<int>> edge_original_nodes = {
		{0, 1, 4},
		{1, 2, 5},
		{2, 3, 6},
		{3, 0, 7}
	};

	for (int edge = 0; edge < 4; edge++) {
		for (int local_node_idx : edge_original_nodes[edge]) {
			int global_node_id = original_nodes[local_node_idx];
			if (original_constraints[global_node_id].empty()) {
				edge_fully_constrained[edge] = false;
				break;
			}
		}
	}

	for (int j = 0; j < nodes_per_side; j++) {
		double eta = gll_points[j];
		for (int i = 0; i < nodes_per_side; i++) {
			double xi = gll_points[i];

			quad8_shape_functions(xi, eta, N);

			double x = 0.0, y = 0.0, z = 0.0;
			for (int k = 0; k < 8; k++) {
				x += N[k] * orig_x[k];
				y += N[k] * orig_y[k];
				z += N[k] * orig_z[k];
			}

			bool node_exists = false;
			int existing_id = -1;

			for (const auto& node : nodes) {
				double dx = node->getX() - x;
				double dy = node->getY() - y;
				double dz = node->getZ() - z;
				double dist_sq = dx * dx + dy * dy + dz * dz;
				if (dist_sq < 1e-8) {
					node_exists = true;
					existing_id = node->getID();
					break;
				}
			}

			if (node_exists) {
				new_node_ids.push_back(existing_id);
			}
			else {
				new_node_ids.push_back(next_node_id);
				std::array<double, 3> coords = { x, y, z };
				auto new_node = std::make_shared<Node>(next_node_id, coords);

				if (is_quad_node_on_constrained_edge(i, j, nodes_per_side, edge_fully_constrained)) {
					apply_constraints_new_quad_nodes(new_node, xi, eta, original_nodes, original_constraints);
				}

				nodes.push_back(new_node);
				next_node_id++;
			}

			new_x.push_back(x);
			new_y.push_back(y);
			new_z.push_back(z);
		}
	}

	elem->set_nodes(new_node_ids);
	elem->set_coords(new_x, new_y, new_z);
}

bool Data::is_quad_node_on_constrained_edge(int i, int j, int nodes_per_side,
	const std::vector<bool>& edge_fully_constrained) {
	if (j == 0 && edge_fully_constrained[0])  
		return true;
	if (i == nodes_per_side - 1 && edge_fully_constrained[1]) 
		return true;
	if (j == nodes_per_side - 1 && edge_fully_constrained[2]) 
		return true;
	if (i == 0 && edge_fully_constrained[3]) 
		return true;

	return false;
}

void Data::apply_constraints_new_quad_nodes(
	std::shared_ptr<Node> new_node, double xi, double eta,
	const std::vector<int>& original_nodes,
	const std::map<int, std::map<int, double>>& original_constraints) {

	int nearest_corner = 0;
	double min_dist = std::numeric_limits<double>::max();

	std::vector<std::pair<double, double>> corner_coords = {
		{-1, -1}, {1, -1}, {1, 1}, {-1, 1}
	};

	for (int corner = 0; corner < 4; corner++) {
		double dist = std::sqrt(std::pow(xi - corner_coords[corner].first, 2) +
			std::pow(eta - corner_coords[corner].second, 2));
		if (dist < min_dist) {
			min_dist = dist;
			nearest_corner = corner;
		}
	}

	int original_node_id = original_nodes[nearest_corner];
	if (original_constraints.count(original_node_id)) {
		const auto& constraints = original_constraints.at(original_node_id);
		for (const auto& constraint : constraints) {
			new_node->set_constraints(constraint.first, constraint.second);
		}
	}
}

void Data::build_spectral_hex_element(std::shared_ptr<Element> elem) {
	int order = elem->get_order();
	int nodes_per_side = order + 1;

	std::vector<int> original_nodes = elem->get_node();
	if (original_nodes.size() != 20) {
		throw std::runtime_error("Expected HEX20 element with 20 nodes for spectral element.");
	}

	std::map<int, std::map<int, double>> original_constraints;
	for (int node_id : original_nodes) {
		auto it = std::find_if(nodes.begin(), nodes.end(),
			[node_id](const std::shared_ptr<Node>& n) { return n->getID() == node_id; });
		if (it != nodes.end()) {
			original_constraints[node_id] = (*it)->constraints;
		}
	}

	std::vector<double> orig_x, orig_y, orig_z;
	for (int node_id : original_nodes) {
		auto it = std::find_if(nodes.begin(), nodes.end(),
			[node_id](const std::shared_ptr<Node>& n) { return n->getID() == node_id; });
		if (it == nodes.end()) {
			throw std::runtime_error("Node " + std::to_string(node_id) + " not found");
		}
		orig_x.push_back((*it)->getX());
		orig_y.push_back((*it)->getY());
		orig_z.push_back((*it)->getZ());
	}

	std::vector<double> gll_points, gll_weights;
	compute_gll_nodes_weights(order, gll_points, gll_weights);

	auto hex20_shape_functions = [](double xi, double eta, double zeta, std::vector<double>& N) {
		N.resize(20);

		N[0] = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta) * (-2 - xi - eta - zeta);
		N[1] = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta) * (-2 + xi - eta - zeta);
		N[2] = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta) * (-2 + xi + eta - zeta);
		N[3] = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta) * (-2 - xi + eta - zeta);
		N[4] = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta) * (-2 - xi - eta + zeta);
		N[5] = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta) * (-2 + xi - eta + zeta);
		N[6] = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta) * (-2 + xi + eta + zeta);
		N[7] = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta) * (-2 - xi + eta + zeta);

		N[8] = 0.25 * (1 - xi * xi) * (1 - eta) * (1 - zeta);
		N[9] = 0.25 * (1 + xi) * (1 - eta * eta) * (1 - zeta);
		N[10] = 0.25 * (1 - xi * xi) * (1 + eta) * (1 - zeta);
		N[11] = 0.25 * (1 - xi) * (1 - eta * eta) * (1 - zeta);
		N[12] = 0.25 * (1 - xi * xi) * (1 - eta) * (1 + zeta);
		N[13] = 0.25 * (1 + xi) * (1 - eta * eta) * (1 + zeta);
		N[14] = 0.25 * (1 - xi * xi) * (1 + eta) * (1 + zeta);
		N[15] = 0.25 * (1 - xi) * (1 - eta * eta) * (1 + zeta);
		N[16] = 0.25 * (1 - xi) * (1 - eta) * (1 - zeta * zeta);
		N[17] = 0.25 * (1 + xi) * (1 - eta) * (1 - zeta * zeta);
		N[18] = 0.25 * (1 + xi) * (1 + eta) * (1 - zeta * zeta);
		N[19] = 0.25 * (1 - xi) * (1 + eta) * (1 - zeta * zeta);
		};

	std::vector<double> N(20);
	std::vector<int> new_node_ids;
	std::vector<double> new_x, new_y, new_z;

	int next_node_id = nodes.size() + 1;
	if (!nodes.empty()) {
		int max_id = 0;
		for (const auto& node : nodes) {
			max_id = std::max(max_id, node->getID());
		}
		next_node_id = max_id + 1;
	}

	std::vector<bool> face_fully_constrained(6, true);
	std::vector<std::vector<int>> face_original_nodes = {
		{0, 1, 2, 3, 8, 9, 10, 11},
		{0, 3, 7, 4, 11, 19, 15, 16},
		{1, 2, 6, 5, 9, 18, 13, 17},
		{2, 3, 7, 6, 10, 19, 14, 18},
		{0, 1, 5, 4, 8, 17, 12, 16},
		{4, 5, 6, 7, 12, 13, 14, 15}
	};
	for (int face = 0; face < 6; face++) {
		for (int local_node_idx : face_original_nodes[face]) {
			int global_node_id = original_nodes[local_node_idx];
			if (original_constraints[global_node_id].empty()) {
				face_fully_constrained[face] = false;
				break;
			}
		}
	}

	for (int k = 0; k < nodes_per_side; k++) {
		double zeta = gll_points[k];
		for (int j = 0; j < nodes_per_side; j++) {
			double eta = gll_points[j];
			for (int i = 0; i < nodes_per_side; i++) {
				double xi = gll_points[i];

				hex20_shape_functions(xi, eta, zeta, N);

				double x = 0.0, y = 0.0, z = 0.0;
				for (int m = 0; m < 20; m++) {
					x += N[m] * orig_x[m];
					y += N[m] * orig_y[m];
					z += N[m] * orig_z[m];
				}

				bool node_exists = false;
				int existing_id = -1;

				for (const auto& node : nodes) {
					double dx = node->getX() - x;
					double dy = node->getY() - y;
					double dz = node->getZ() - z;
					double dist_sq = dx * dx + dy * dy + dz * dz;
					if (dist_sq < 1e-8) {
						node_exists = true;
						existing_id = node->getID();
						break;
					}
				}

				if (node_exists) {
					new_node_ids.push_back(existing_id);
				}
				else {
					new_node_ids.push_back(next_node_id);
					std::array<double, 3> coords = { x, y, z };
					auto new_node = std::make_shared<Node>(next_node_id, coords);

					if (is_hex_node_on_constrained_face(i, j, k, nodes_per_side, face_fully_constrained)) {
						apply_constraints_new_hex_nodes(new_node, xi, eta, zeta, original_nodes, original_constraints);
					}

					nodes.push_back(new_node);
					next_node_id++;
				}

				new_x.push_back(x);
				new_y.push_back(y);
				new_z.push_back(z);
			}
		}
	}

	elem->set_nodes(new_node_ids);
	elem->set_coords(new_x, new_y, new_z);
}

bool Data::is_hex_node_on_constrained_face(
	int i, int j, int k, int nodes_per_side,
	const std::vector<bool>& face_fully_constrained) {
	if (j == 0 && face_fully_constrained[0])
		return true;
	if (i == 0 && face_fully_constrained[1])
		return true;
	if (i == nodes_per_side - 1 && face_fully_constrained[2])
		return true;
	if (j == nodes_per_side - 1 && face_fully_constrained[3])
		return true;
	if (k == nodes_per_side - 1 == 0 && face_fully_constrained[4])
		return true;
	if (k == 0 && face_fully_constrained[5])
		return true;

	return false;
}

void Data::apply_constraints_new_hex_nodes(
	std::shared_ptr<Node> new_node,
	double xi, double eta, double zeta,
	const std::vector<int>& original_nodes,
	const std::map<int, std::map<int, double>>& original_constraints) {
	int nearest_corner = 0;
	double min_dist = std::numeric_limits<double>::max();

	std::vector<std::tuple<double, double, double>> corner_coords = {
		{-1, -1, -1},
		{ 1, -1, -1},
		{ 1,  1, -1},
		{-1,  1, -1},
		{-1, -1,  1},
		{ 1, -1,  1},
		{ 1,  1,  1},
		{-1,  1,  1}
	};

	for (int corner = 0; corner < 8; corner++) {
		auto [cx, cy, cz] = corner_coords[corner];
		double dist = std::sqrt(std::pow(xi - cx, 2) +
			std::pow(eta - cy, 2) +
			std::pow(zeta - cz, 2));
		if (dist < min_dist) {
			min_dist = dist;
			nearest_corner = corner;
		}
	}

	int original_node_id = original_nodes[nearest_corner];
	if (original_constraints.count(original_node_id)) {
		const auto& constraints = original_constraints.at(original_node_id);
		for (const auto& [dof, value] : constraints) {
			new_node->set_constraints(dof, value);
		}
	}
}

void Data::renumber_nodes() {
	logger& log = logger::log();
	log.print("Starting node renumbering");

	std::set<int> used_nodes;
	for (const auto& elem : elements) {
		const auto& elem_nodes = elem->get_node();
		used_nodes.insert(elem_nodes.begin(), elem_nodes.end());
	}

	std::map<int, int> old_to_new;
	int new_id = 1;
	for (int old_id : used_nodes) {
		old_to_new[old_id] = new_id++;
	}

	for (auto& elem : elements) {
		std::vector<int> new_nodes;
		for (int old_node : elem->get_node()) {
			new_nodes.push_back(old_to_new.at(old_node));
		}
		elem->set_nodes(new_nodes);
	}

	std::vector<std::shared_ptr<Node>> new_nodes_list;
	new_nodes_list.resize(old_to_new.size());

	for (const auto& old_node : nodes) {
		int old_id = old_node->getID();
		if (old_to_new.count(old_id)) {
			int new_id = old_to_new[old_id];
			auto new_node = std::make_shared<Node>(new_id,
				std::array<double, 3>{ old_node->getX(), old_node->getY(), old_node->getZ()	}
			);
			new_node->load = old_node->load;
			new_node->constraints = old_node->constraints;
			new_node->results = old_node->results;

			new_nodes_list[new_id - 1] = new_node;
		}
	}

	nodes = std::move(new_nodes_list);

	log.print("Node renumbering completed. Total nodes: " + std::to_string(nodes.size()));
}