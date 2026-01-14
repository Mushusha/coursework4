#include "Data.h"

bool Data::is_dynamic = false;
double Data::reference_element_size = 1.0;
double Data::damping = 0.0;

Data::Data(const Data& other)
	: dim(other.dim),
	filename(other.filename),
	analisys_type(other.analisys_type),
	max_time(other.max_time),
	max_iter(other.max_iter),
	iter_res_output(other.iter_res_output),
	omega(other.omega),
	Amp(other.Amp),
	elements(other.elements),
	nodes(other.nodes),
	node_id_old_to_new(other.node_id_old_to_new) {
}
Data& Data::operator=(const Data& other) {
	if (this != &other) {
		dim = other.dim;
		filename = other.filename;
		analisys_type = other.analisys_type;
		max_time = other.max_time;
		max_iter = other.max_iter;
		iter_res_output = other.iter_res_output;
		omega = other.omega;
		Amp = other.Amp;
		elements = other.elements;
		nodes = other.nodes;
		node_id_old_to_new = other.node_id_old_to_new;
	}
	return *this;
}
Data::Data(Data&& other) noexcept
	: dim(std::exchange(other.dim, 0)),
	filename(std::move(other.filename)),
	analisys_type(std::exchange(other.analisys_type, 0)),
	max_time(std::exchange(other.max_time, 0)),
	max_iter(std::exchange(other.max_iter, 0)),
	iter_res_output(std::exchange(other.iter_res_output, 0)),
	omega(std::exchange(other.omega, 0)),
	Amp(std::exchange(other.Amp, 0)),
	elements(std::move(other.elements)),
	nodes(std::move(other.nodes)),
	node_id_old_to_new(std::move(other.node_id_old_to_new)) {
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
		node_id_old_to_new = std::move(other.node_id_old_to_new);
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
	max_time = parser->settings.max_time;
	max_iter = parser->settings.max_iter;
	iter_res_output = parser->settings.iter_res_output;

	omega = parser->settings.omega; // this params identical for all components
	Amp = parser->settings.Amp;
	
	is_dynamic = parser->settings.analisys_type == "dynamic";
	damping = parser->settings.d;
	
	create_nodes(parser);
	create_constraints(parser);
	create_elements(parser);
	create_constants(parser);
	create_load(parser);
	create_D(parser);
	
	compute_reference_element_size();

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
		else if (type == HEXSEM && order > 2) {
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
			case PYR:
				elem = std::make_shared<Pyr>(Pyr(parser->mesh.elem_id[i], PYR, elem_nodes));
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

	std::map<int, int> node_id_to_index;
	for (size_t idx = 0; idx < nodes.size(); ++idx) {
		node_id_to_index[nodes[idx]->getID()] = static_cast<int>(idx);
	}

	auto get_node_by_id = [&](int node_id) -> std::shared_ptr<Node> {
		auto it = node_id_to_index.find(node_id);
		if (it == node_id_to_index.end()) {
			throw std::runtime_error("Node ID " + std::to_string(node_id) + " not found in node_id_to_index map");
		}
		return nodes[it->second];
	};

	auto set_pole = [](auto* elem, const std::array<double, 3>& point) {
		elem->pole_x = point[0];
		elem->pole_y = point[1];
		elem->pole_z = point[2];
	};

	for (const auto& inf : parser->infinite) {
		for (int i = 0; i < inf.size; i++) {
			const int elem_id = inf.apply_to[2 * i];
			const int side = inf.apply_to[2 * i + 1];
			
			auto parent_elem = this->elements[elem_id - 1];
			ElemType parent_type = parent_elem->get_type();
			int parent_order = parent_elem->get_order();
			int NODES = parent_order + 1;

			std::vector<int> side_nodes = parent_elem->edge_to_node(side);
			std::vector<int> boundary_glob_nodes;
			for (int& local_node : side_nodes) {
				boundary_glob_nodes.push_back(parent_elem->get_node(local_node));
			}

			std::vector<double> gr_x, gr_w;
			compute_gauss_radau_nodes_weights(NODES, gr_x, gr_w);

			std::vector<std::vector<int>> inf_layer_nodes(NODES);
			inf_layer_nodes[0] = boundary_glob_nodes;

			for (int k = 1; k < NODES; ++k) {
				for (size_t n = 0; n < boundary_glob_nodes.size(); ++n) {
					int boundary_node_id = boundary_glob_nodes[n];
					auto boundary_node = get_node_by_id(boundary_node_id);
					
					double t = (gr_x[k] + 1.0) / 2.0;
					std::array<double, 3> coords;
					for (int j = 0; j < 3; j++) {
						double x_boundary = boundary_node->getCoord(j);
						double x_pole = inf.point[j];
						double x_inf = 2.0 * x_boundary - x_pole;
						coords[j] = x_boundary + (x_inf - x_boundary) * t;
					}

					int max_node_id = 0;
					for (const auto& node : nodes) {
						max_node_id = std::max(max_node_id, node->getID());
					}
					int new_node_id = max_node_id + 1;
					
					this->nodes.push_back(std::make_shared<Node>(new_node_id, coords));
					node_id_to_index[new_node_id] = static_cast<int>(nodes.size() - 1);
					inf_layer_nodes[k].push_back(new_node_id);
				}
			}

			std::vector<int> inf_elem_nodes;
			int num_nodes_per_layer = static_cast<int>(boundary_glob_nodes.size());
			for (int j = 0; j < num_nodes_per_layer; ++j) {
				for (int i = 0; i < NODES; ++i) {
					if (i < static_cast<int>(inf_layer_nodes.size()) && 
					    j < static_cast<int>(inf_layer_nodes[i].size())) {
						inf_elem_nodes.push_back(inf_layer_nodes[i][j]);
					}
				}
			}

			std::vector<double> x1, y1, z1;
			for (auto& node_id : inf_elem_nodes) {
				auto node = get_node_by_id(node_id);
				x1.push_back(node->getX());
				y1.push_back(node->getY());
				z1.push_back(node->getZ());
			}

			std::shared_ptr<Element> new_elem;
			
			if (parent_type == QUADSEM || parent_type == HEXSEM) {
				#define CREATE_SPECTRAL_INF_ELEM(ElemType, ElemTypeEnum, N) \
					new_elem = std::make_shared<ElemType<N>>( \
						ElemType<N>(this->elements.size() + 1, ElemTypeEnum, inf_elem_nodes)); \
					new_elem->set_coords(x1, y1, z1); \
					new_elem->set_order(parent_order); \
					auto* elem_ptr = dynamic_cast<ElemType<N>*>(new_elem.get()); \
					set_pole(elem_ptr, inf.point); \
					if (is_dynamic) { \
						auto* inf_elem = dynamic_cast<ElemType<N>*>(new_elem.get()); \
						inf_elem->omega = omega; \
						inf_elem->is_dynamic = true; \
					}

				if (dim == 2) {
					switch (NODES) {
						case 2: { CREATE_SPECTRAL_INF_ELEM(SpectralInfQuad, INFQUADSEM, 2); break; }
						case 3: { CREATE_SPECTRAL_INF_ELEM(SpectralInfQuad, INFQUADSEM, 3); break; }
						case 4: { CREATE_SPECTRAL_INF_ELEM(SpectralInfQuad, INFQUADSEM, 4); break; }
						case 5: { CREATE_SPECTRAL_INF_ELEM(SpectralInfQuad, INFQUADSEM, 5); break; }
						case 6: { CREATE_SPECTRAL_INF_ELEM(SpectralInfQuad, INFQUADSEM, 6); break; }
						default: throw std::runtime_error("Unsupported order for SpectralInfQuad");
					}
				}
				else {
					switch (NODES) {
						case 2: { CREATE_SPECTRAL_INF_ELEM(SpectralInfHex, INFHEXSEM, 2); break; }
						case 3: { CREATE_SPECTRAL_INF_ELEM(SpectralInfHex, INFHEXSEM, 3); break; }
						case 4: { CREATE_SPECTRAL_INF_ELEM(SpectralInfHex, INFHEXSEM, 4); break; }
						case 5: { CREATE_SPECTRAL_INF_ELEM(SpectralInfHex, INFHEXSEM, 5); break; }
						case 6: { CREATE_SPECTRAL_INF_ELEM(SpectralInfHex, INFHEXSEM, 6); break; }
						default: throw std::runtime_error("Unsupported order for SpectralInfHex");
					}
				}
				#undef CREATE_SPECTRAL_INF_ELEM
			}
			else if (dim == 2) {
				std::vector<int> inf_quad_nodes = { 
					inf_layer_nodes[0][0], inf_layer_nodes[1][0], 
					inf_layer_nodes[1][1], inf_layer_nodes[0][1] 
				};

				std::vector<double> xq, yq, zq;
				for (auto& node_id : inf_quad_nodes) {
					auto node = get_node_by_id(node_id);
					xq.push_back(node->getX());
					yq.push_back(node->getY());
					zq.push_back(node->getZ());
				}

				new_elem = std::make_shared<InfQuad>(InfQuad(this->elements.size() + 1, INFQUAD, inf_quad_nodes));
				new_elem->set_coords(xq, yq, zq);
				new_elem->set_order(parent_order);
				
				auto* iq = dynamic_cast<InfQuad*>(new_elem.get());
				set_pole(iq, inf.point);
				
				if (is_dynamic) {
					iq->omega = omega;
				}
			}
			else {
				std::vector<int> inf_hex_nodes;
				for (auto& node : inf_layer_nodes[0])
					inf_hex_nodes.push_back(node);
				for (auto& node : inf_layer_nodes[1])
					inf_hex_nodes.push_back(node);

				std::vector<double> xh, yh, zh;
				for (auto& node_id : inf_hex_nodes) {
					auto node = get_node_by_id(node_id);
					xh.push_back(node->getX());
					yh.push_back(node->getY());
					zh.push_back(node->getZ());
				}

				new_elem = std::make_shared<InfHex>(InfHex(this->elements.size() + 1, INFHEX, inf_hex_nodes));
				new_elem->set_coords(xh, yh, zh);
				new_elem->set_order(parent_order);

				auto* ih = dynamic_cast<InfHex*>(new_elem.get());
				set_pole(ih, inf.point);

				if (is_dynamic) {
					ih->omega = omega;
				}
			}

			transfer_pressure_to_infinite_element(parser, parent_elem, new_elem, elem_id, side);

			this->elements.push_back(new_elem);
			num_inf_elems++;
		}
	}
}

void Data::transfer_pressure_to_infinite_element(std::shared_ptr<const Parser> parser,
	std::shared_ptr<Element> parent_elem, std::shared_ptr<Element> inf_elem,
	int parent_elem_id, int parent_side) {
	
	double pressure_value = 0.0;
	int pressure_type = PRESSURE;
	bool found_pressure = false;
	
	for (const auto& load_entry : parser->load) {
		if (load_entry.inf) continue;
		if (load_entry.type != PRESSURE && load_entry.type != PRESSURESEM) continue;
		
		for (size_t j = 0; j < load_entry.apply_to.size() / 2; j++) {
			int load_elem_id = load_entry.apply_to[2 * j];
			int load_side = load_entry.apply_to[2 * j + 1];
			
			if (load_elem_id == parent_elem_id && load_side == parent_side) {
				pressure_value = load_entry.data[0];
				pressure_type = load_entry.type;
				found_pressure = true;
				break;
			}
		}
		if (found_pressure) break;
	}
	
	if (!found_pressure) {
		return;
	}
	
	int boundary_edge = -1;
	if (this->dim == 2) {
		boundary_edge = 3;
	} else {
		boundary_edge = 4;
	}
	
	std::vector<int> parent_edge_nodes = parent_elem->edge_to_node(parent_side);
	std::vector<int> inf_edge_nodes;
	
	const double tol = 1e-6;
	
	if (inf_elem->get_type() == INFQUADSEM || (inf_elem->get_type() == INFQUAD && parent_edge_nodes.size() > 2)) {
		inf_edge_nodes.clear();
		for (size_t i = 0; i < parent_edge_nodes.size(); i++) {
			double px = parent_elem->get_coord(parent_edge_nodes[i], 0);
			double py = parent_elem->get_coord(parent_edge_nodes[i], 1);
			double pz = parent_elem->get_coord(parent_edge_nodes[i], 2);
			
			for (int j = 0; j < inf_elem->nodes_count(); j++) {
				double ix = inf_elem->get_coord(j, 0);
				double iy = inf_elem->get_coord(j, 1);
				double iz = inf_elem->get_coord(j, 2);
				
				double dx = std::abs(px - ix);
				double dy = std::abs(py - iy);
				double dz = std::abs(pz - iz);
				
				if (dx < tol && dy < tol && dz < tol) {
					if (std::find(inf_edge_nodes.begin(), inf_edge_nodes.end(), j) == inf_edge_nodes.end()) {
						inf_edge_nodes.push_back(j);
						break;
					}
				}
			}
		}
	} else {
		inf_edge_nodes = inf_elem->edge_to_node(boundary_edge);
	}
	
	parent_elem->remove_load_on_edge(parent_side);
	
	if (boundary_edge >= 0) {
		inf_elem->set_load(pressure_type, boundary_edge, {pressure_value, 0.0, 0.0, 0.0, 0.0, 0.0});
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
						for (int node = 0; node < load.apply_to.size(); node += 2) {
							int old_node_id = load.apply_to[node];
							int new_node_id = node_id_old_to_new.count(old_node_id) 
								? node_id_old_to_new[old_node_id] 
								: old_node_id;
							nodes[new_node_id - 1]->load[i] = load.data[i];
						}
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

	int total_nodes = nodes_per_side * nodes_per_side;
	new_node_ids.reserve(total_nodes);
	new_x.reserve(total_nodes);
	new_y.reserve(total_nodes);
	new_z.reserve(total_nodes);
	
	struct NodeData {
		double x, y, z;
		double xi, eta;
		int i, j;
		bool needs_constraints;
		int existing_id;
		bool node_exists;
	};
	
	std::vector<NodeData> node_data_list(total_nodes);
	
	size_t initial_nodes_size = nodes.size();
	
	std::vector<int> indices(total_nodes);
	std::iota(indices.begin(), indices.end(), 0);
	
	std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int idx) {
		int j = idx / nodes_per_side;
		int i = idx % nodes_per_side;
		double eta = gll_points[j];
		double xi = gll_points[i];
		
		std::vector<double> N_local(8);
		quad8_shape_functions(xi, eta, N_local);
		
		double x = 0.0, y = 0.0, z = 0.0;
		for (int k = 0; k < 8; k++) {
			x += N_local[k] * orig_x[k];
			y += N_local[k] * orig_y[k];
			z += N_local[k] * orig_z[k];
		}
		
		bool node_exists = false;
		int existing_id = -1;
		const double tol_sq = 1e-8;
		
		for (size_t n_idx = 0; n_idx < initial_nodes_size; ++n_idx) {
			const auto& node = nodes[n_idx];
			double dx = node->getX() - x;
			double dy = node->getY() - y;
			double dz = node->getZ() - z;
			double dist_sq = dx * dx + dy * dy + dz * dz;
			if (dist_sq < tol_sq) {
				node_exists = true;
				existing_id = node->getID();
				break;
			}
		}
		
		node_data_list[idx].x = x;
		node_data_list[idx].y = y;
		node_data_list[idx].z = z;
		node_data_list[idx].xi = xi;
		node_data_list[idx].eta = eta;
		node_data_list[idx].i = i;
		node_data_list[idx].j = j;
		node_data_list[idx].needs_constraints = is_quad_node_on_constrained_edge(i, j, nodes_per_side, edge_fully_constrained);
		node_data_list[idx].existing_id = existing_id;
		node_data_list[idx].node_exists = node_exists;
	});
	
	for (const auto& data : node_data_list) {
		new_x.push_back(data.x);
		new_y.push_back(data.y);
		new_z.push_back(data.z);
		
		if (data.node_exists) {
			new_node_ids.push_back(data.existing_id);
		}
		else {
			new_node_ids.push_back(next_node_id);
			std::array<double, 3> coords = { data.x, data.y, data.z };
			auto new_node = std::make_shared<Node>(next_node_id, coords);
			
			if (data.needs_constraints) {
				apply_constraints_new_quad_nodes(new_node, data.xi, data.eta, original_nodes, original_constraints);
			}
			
			nodes.push_back(new_node);
			next_node_id++;
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

	int total_nodes = nodes_per_side * nodes_per_side * nodes_per_side;
	new_node_ids.reserve(total_nodes);
	new_x.reserve(total_nodes);
	new_y.reserve(total_nodes);
	new_z.reserve(total_nodes);
	
	struct HexNodeData {
		double x, y, z;
		double xi, eta, zeta;
		int i, j, k;
		bool needs_constraints;
		int existing_id;
		bool node_exists;
	};
	
	std::vector<HexNodeData> node_data_list(total_nodes);
	
	size_t initial_nodes_size = nodes.size();
	
	std::vector<int> indices(total_nodes);
	std::iota(indices.begin(), indices.end(), 0);
	
	std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int idx) {
		int k = idx / (nodes_per_side * nodes_per_side);
		int remainder = idx % (nodes_per_side * nodes_per_side);
		int j = remainder / nodes_per_side;
		int i = remainder % nodes_per_side;
		
		double zeta = gll_points[k];
		double eta = gll_points[j];
		double xi = gll_points[i];
		
		std::vector<double> N_local(20);
		hex20_shape_functions(xi, eta, zeta, N_local);
		
		double x = 0.0, y = 0.0, z = 0.0;
		for (int m = 0; m < 20; m++) {
			x += N_local[m] * orig_x[m];
			y += N_local[m] * orig_y[m];
			z += N_local[m] * orig_z[m];
		}
		
		bool node_exists = false;
		int existing_id = -1;
		const double tol_sq = 1e-8;
		
		for (size_t n_idx = 0; n_idx < initial_nodes_size; ++n_idx) {
			const auto& node = nodes[n_idx];
			double dx = node->getX() - x;
			double dy = node->getY() - y;
			double dz = node->getZ() - z;
			double dist_sq = dx * dx + dy * dy + dz * dz;
			if (dist_sq < tol_sq) {
				node_exists = true;
				existing_id = node->getID();
				break;
			}
		}
		
		node_data_list[idx].x = x;
		node_data_list[idx].y = y;
		node_data_list[idx].z = z;
		node_data_list[idx].xi = xi;
		node_data_list[idx].eta = eta;
		node_data_list[idx].zeta = zeta;
		node_data_list[idx].i = i;
		node_data_list[idx].j = j;
		node_data_list[idx].k = k;
		node_data_list[idx].needs_constraints = is_hex_node_on_constrained_face(i, j, k, nodes_per_side, face_fully_constrained);
		node_data_list[idx].existing_id = existing_id;
		node_data_list[idx].node_exists = node_exists;
	});
	
	for (const auto& data : node_data_list) {
		new_x.push_back(data.x);
		new_y.push_back(data.y);
		new_z.push_back(data.z);
		
		if (data.node_exists) {
			new_node_ids.push_back(data.existing_id);
		}
		else {
			new_node_ids.push_back(next_node_id);
			std::array<double, 3> coords = { data.x, data.y, data.z };
			auto new_node = std::make_shared<Node>(next_node_id, coords);
			
			if (data.needs_constraints) {
				apply_constraints_new_hex_nodes(new_node, data.xi, data.eta, data.zeta, original_nodes, original_constraints);
			}
			
			nodes.push_back(new_node);
			next_node_id++;
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
	if (k == nodes_per_side - 1 && face_fully_constrained[4])
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

	node_id_old_to_new = old_to_new;

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

void Data::compute_reference_element_size() {
	double total_size = 0.0;
	int finite_elem_count = 0;
	
	for (int i = 0; i < elements_count() - num_inf_elems; i++) {
		double vol = elements[i]->Volume();
		if (vol > 0.0) {
			if (dim == 2) {
				total_size += std::sqrt(vol);
			} else {
				total_size += std::cbrt(vol);
			}
			finite_elem_count++;
		}
	}
	
	if (finite_elem_count > 0) {
		reference_element_size = total_size / finite_elem_count;
	} else {
		reference_element_size = 1.0;
	}
	
	if (reference_element_size < 1e-10) reference_element_size = 1e-10;
	if (reference_element_size > 1e10) reference_element_size = 1e10;
}