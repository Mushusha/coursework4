#include "Data.h"

Data::Data(std::shared_ptr <Parser> p) : parser(p) {
	this->dim = 2; // need read from fc
	create_nodes();
	create_elements();
	create_constrains();
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
}

void Data::create_constrains() {
	for (int id = 0; id < parser->restraints.size(); id++)
		for (int node = 0; node < parser->restraints[id].size; node++)
			for (int i = 0; i < 6; i++)
				if (parser->restraints[id].flag[i])
					nodes[node].set_constrains(i, parser->restraints[id].data[i]);
}

void Data::fillGlobalK() {
	int inf_count = 0; // ??
	int nodes_count = parser->mesh.nodes_count;
	int elems_count = parser->mesh.elems_count;

	Eigen::SparseMatrix <double> K(dim * (nodes_count + inf_count), dim * (nodes_count + inf_count));
	std::vector <Eigen::Triplet <double>> tripl_vec;
	for (int i = 0; i < elems_count; i++)
		for (int j = 0; j < elements[i]->nodes_count() * dim; j++)
			for (int k = 0; k < elements[i]->nodes_count() * dim; k++) {
				Eigen::MatrixXd loc_k = elements[i]->locK();
				Eigen::Triplet <double> trpl(dim * (elements[i]->get_nodes(j / dim) - 1) + j % dim, dim * (elements[i]->get_nodes(k / dim) - 1) + k % dim, loc_k(j, k));
				tripl_vec.push_back(trpl);
			}

	K.setFromTriplets(tripl_vec.begin(), tripl_vec.end());	
} 

void Data::fillGlobalLoad() {

}

void Data::fillConstrains() {
	for (int node = 0; node < nodes.size(); node++)
		for (auto const& c: nodes[node].constrains)
			for (int i = 0; i < K.outerSize(); i++) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
					if ((it.row() == (node - 1) * dim + c.first) ||
						(it.col() == (node - 1) * dim + c.first) && (it.row() != it.col()))
						it.valueRef() = 0.0;
				R.coeffRef(i) *= c.second;
			}
}

void Data::addToGlobalK(int first_index, int second_index, double value) {
	Eigen::Triplet <double> tripl (first_index, second_index, value);
	K.setFromTriplets(&tripl, &tripl + 1);
}

void Data::addToGlobalR() {
}
