#include "Element.h"


double Young = 2e+11;
double Poisson = 0.3; // need read from fc


void Element::set_coords(std::vector <double> x, std::vector <double> y, std::vector <double> z) {
	for (int i = 0; i < x.size(); i++) {
		this->x.push_back(x[i]);
		this->y.push_back(y[i]);
		this->z.push_back(z[i]);
	}
}

void Element::set_load(int edge, int comp, double value) {
	std::pair <int, int> pair(edge, comp);
	load.insert({ pair, value });
}

Eigen::MatrixXd Element::twoMatrixD() {

	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 3);
	D << 1, Poisson, 0,
		Poisson, 1, 0,
		0, 0, (1 - Poisson) / 2;

	D *= Young / (1 - pow(Poisson, 2));

	return D;
}

Eigen::MatrixXd Element::threeMatrixD() {
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(6, 6);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 3; j++) {
			if (i == j && i < 3)
				D(i, j) = 1;
			if (i != j && i < 3)
				D(i, j) = Poisson / (1 - Poisson);
			if (i >= 3)
				D(i, i) = (1 - 2 * Poisson) / (2 * (1 - Poisson));
		}

	D *= Young * (1 - Poisson) / ((1 + Poisson) * (1 - 2 * Poisson));
	return D;
}



//Eigen::SparseMatrix <double> globalK(std::shared_ptr <Parser> p, int numNodes, int dim) {
//	int infCount = 0; //
//	int nodesCount = p->mesh.nodes_count;
//	int numElements = p->mesh.elems_count;
//	std::vector <std::shared_ptr <Element>> elem;
//	elem = createElements(p, numNodes);
//	Eigen::SparseMatrix <double> K(dim * (nodesCount + infCount), dim * (nodesCount + infCount));
//	std::vector <Eigen::Triplet <double>> tripl;
//	for (int i = 0; i < numElements; i++)
//		for (int j = 0; j < numNodes * dim; j++)
//			for (int k = 0; k < numNodes * dim; k++) {
//				Eigen::Triplet <double> trpl(dim * (elem[i]->nodes[j / dim] - 1) + j % dim, dim * (elem[i]->nodes[k / dim] - 1) + k % dim, elem[i]->locK(j, k));
//				tripl.push_back(trpl);
//			}
//
//	K.setFromTriplets(tripl.begin(), tripl.end());
//
//	for (int i = 0; i < K.outerSize(); ++i)
//		for (Eigen::SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
//			for (int j = 0; j < p->restraints.size(); j++)
//				for (int k = 0; k < p->restraints[j].apply_to.size(); k++) 
//					if (it.row() == indexRestrain(p, j, k, dim) || it.col() == indexRestrain(p, j, k, dim)) 
//						it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
//						
//	K.makeCompressed();
//	return K;
//}
//
//template <const int NODES, const int DIM>
//Eigen::SparseVector <double> globalLoad(std::shared_ptr <Parser> p, int numNodes, int dim) {
//	int infCount = 0; //
//	int nodesCount = p->mesh.nodes_count;
//	int numElements = p->mesh.elems_count;
//	std::vector <std::shared_ptr <Element>> elem;
//
//	Eigen::SparseVector <double> L(dim * (nodesCount + infCount));
//	for (int i = 0; i < numElements; i++)
//		for (int j = 0; j < numNodes * dim; j++)
//			L.insert(dim * (elem[i]->nodes[j / dim] - 1) + j % dim) = 1; // тут лоады дописать
//	return L;
//}
//
//int indexRestrain(std::shared_ptr <Parser> p, int count, int apply, int dim) {
//	int temp = 0;
//	for (int i = 0; i < 6; i++)
//		if (p->restraints[count].flag[i] == 1) {
//			temp = i;
//			break;
//		}
//	return dim * (p->restraints[count].apply_to[apply] - 1) + temp;
//}
//
//template <const int NODES, const int DIM>
//std::vector <std::shared_ptr <Element <NODES, DIM>>> createElements(std::shared_ptr <Parser> p, int numNodes) {
//	int numElements = p->mesh.elems_count;
//	std::vector <std::shared_ptr <Element>> elem;
//	for (int i = 0; i < numElements; i++) {
//		std::vector <int> nodes;
//		std::vector <double> x;
//		std::vector <double> y;
//		std::vector <double> z;
//		for (int j = 0; j < numNodes; j++) {
//			nodes.push_back(p->mesh.elem_nodes[numNodes * i + j]);
//			x.push_back(p->mesh.nodes_coord[3 * (nodes[j] - 1)]);
//			y.push_back(p->mesh.nodes_coord[3 * (nodes[j] - 1) + 1]);
//			z.push_back(p->mesh.nodes_coord[3 * (nodes[j] - 1) + 2]);
//		}
//		switch (numNodes) {
//		case 4:
//			elem.push_back(std::make_shared <quadElement>(nodes, x, y));
//			break;
//		case 3:
//			elem.push_back(std::make_shared <triElement>(nodes, x, y));
//			break;
//		default:
//			std::cout << "Error" << endl;
//			break;
//		}
//	}
//	return elem;
//}

// createInf ()