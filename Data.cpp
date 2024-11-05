#include "Data.h"

Data::Data(std::shared_ptr <Parser> p) : parser(p) {
	this->dim = 2; // need read from fc
	// create nodes
	for (int i = 0; i < p->mesh.nodes_count; i++) {
		std::array <double, 3> coords;
		for (int j = 0; j < 3; j++)
			coords[j] = p->mesh.nodes_coord[3 * i + j];
		this->nodes.push_back(Node(p->mesh.node_id[i], coords));
	}
	// create elements
	int node_tmp = 0;
	for (int i = 0; i < p->mesh.elems_count; i++) {
		std::vector<int> elem_nodes;
		std::shared_ptr<Element> elem;
		switch (p->mesh.elem_types[i]) {
		case ElemType::TRI:
			elem_nodes.resize(3);
			for (int j = 0; j < 3; j++)
				elem_nodes[j] = p->mesh.elem_nodes[node_tmp + j];
			node_tmp += 3;
			elem = std::make_shared<triElement>(triElement(p->mesh.elem_id[i], ElemType::TRI, elem_nodes));
			break;
		case ElemType::QUAD:
			elem_nodes.resize(4);
			for (int j = 0; j < 4; j++)
				elem_nodes[j] = p->mesh.elem_nodes[node_tmp + j];
			node_tmp += 4;
			elem = std::make_shared<quadElement>(quadElement(p->mesh.elem_id[i], ElemType::QUAD, elem_nodes));
			break;
		case ElemType::TETRA:
			elem_nodes.resize(4);
			for (int j = 0; j < 4; j++)
				elem_nodes[j] = p->mesh.elem_nodes[node_tmp + j];
			node_tmp += 4;
			elem = std::make_shared<tetraElement>(tetraElement(p->mesh.elem_id[i], ElemType::TETRA, elem_nodes));
			break;
		case ElemType::HEX:
			elem_nodes.resize(8);
			for (int j = 0; j < 8; j++)
				elem_nodes[j] = p->mesh.elem_nodes[node_tmp + j];
			node_tmp += 8;
			elem = std::make_shared<hexElement>(hexElement(p->mesh.elem_id[i], ElemType::HEX, elem_nodes));
			break;
		default:
			std::cout << "Error: incorrect element type";  //need throw
			break;
		}
		std::vector<double> x, y, z;
		for (int i = 0; i < elem_nodes.size(); i++) {
			x.push_back(p->mesh.nodes_coord[3 * elem_nodes[i]]);
			y.push_back(p->mesh.nodes_coord[3 * elem_nodes[i] + 1]);
			z.push_back(p->mesh.nodes_coord[3 * elem_nodes[i] + 2]); // may be Nodes 
		}
		elem->set_coords(x, y, z);
		this->elements.push_back(elem);
	}
}

void Data::fillGlobalK() {
	int infCount = 0; //
	int nodes_count = parser->mesh.nodes_count;
	int elems_count = parser->mesh.elems_count;

	Eigen::SparseMatrix <double> K(dim * (nodes_count + infCount), dim * (nodes_count + infCount));
	std::vector <Eigen::Triplet <double>> tripl_vec;
	for (int i = 0; i < elems_count; i++)
		for (int j = 0; j < elements[i]->nodes_count() * dim; j++)
			for (int k = 0; k < elements[i]->nodes_count() * dim; k++) {
				Eigen::MatrixXd loc_k = elements[i]->locK();
				Eigen::Triplet <double> trpl(dim * (elements[i]->get_nodes(j / dim) - 1) + j % dim, dim * (elements[i]->get_nodes(k / dim) - 1) + k % dim, loc_k(j, k));
				tripl_vec.push_back(trpl);
			}

	K.setFromTriplets(tripl_vec.begin(), tripl_vec.end());

	//for (int i = 0; i < K.outerSize(); ++i)
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
	//		for (int j = 0; j < p->restraints.size(); j++)
	//			for (int k = 0; k < p->restraints[j].apply_to.size(); k++) 
	//				if (it.row() == indexRestrain(p, j, k, dim) || it.col() == indexRestrain(p, j, k, dim)) 
	//					it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
						
 }

void Data::fillGlobalLoad() {
}

void Data::addConstrains() {
}

void Data::addToGlobalK(int first_index, int second_index, double value) {
	Eigen::Triplet <double> tripl (first_index, second_index, value);
	K.setFromTriplets(&tripl, &tripl + 1);
}

void Data::addToGlobalR() {
}
