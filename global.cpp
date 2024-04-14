#include "global.h"
#include "QuadElements.h"
#include "TriangleElements.h"
#include "TetrahedronElements.h"
#include "HexahedronElements.h"

float Young = 20000;
float Poisson = 0.3;
int numNodes;

// Matrix D

Eigen::Matrix3f twoMatrixD() {

	Eigen::Matrix3f D;
	D << 1, Poisson, 0,
		Poisson, 1, 0,
		0, 0, (1 - Poisson) / 2;

	D *= Young / (1 - pow(Poisson, 2));

	return D;
}

Eigen::MatrixXf threeMatrixD() {
	Eigen::MatrixXf D = Eigen::MatrixXf::Zero(6, 6);
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

// globalK

Eigen::SparseMatrix <float> globalK(std::shared_ptr <Parser> p, int numNodes, int dim) {
	int numInf = 0; //
	int numElements = p->mesh.elems_count;
	createElements(p, numNodes);
	Eigen::SparseMatrix <float> K(dim * (numElements + numInf), dim * (numElements + numInf));
	std::vector <Eigen::Triplet <float>> tripl;

	return K;
}

// create

std::vector <std::shared_ptr <Element>> createElements(std::shared_ptr <Parser> p, int numNodes) {
	int numElements = p->mesh.elems_count;
	std::vector <std::shared_ptr <Element>> elem;
	for (int i = 0; i < numElements; i++) {
		std::vector <int> nodes;
		std::vector <float> x;
		std::vector <float> y;
		std::vector <float> z;
		for (int j = 0; j < numNodes; j++) {
			nodes.push_back(p->mesh.elem_nodes[numNodes * i + j]);
			x.push_back(p->mesh.nodes_coord[3 * (nodes[j] - 1)]);
			y.push_back(p->mesh.nodes_coord[3 * (nodes[j] - 1) + 1]);
			z.push_back(p->mesh.nodes_coord[3 * (nodes[j] - 1) + 2]);
		}
		switch (numNodes) {
		case 4: 
			elem.push_back(std::make_shared <quadElement>(nodes, x, y));
			break;
		case 3: 
			elem.push_back(std::make_shared <triangleElement>(nodes, x, y));
			break;
		default: 
			std::cout << "Error" << endl;
			break;
		}
	}
	return elem;
}

// createInf ()