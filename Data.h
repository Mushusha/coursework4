#pragma once

#include "Parser.h"
#include "log.h"
#include "Enums.h"
#include "Node.h"
#include "Element.h"
#include "TriangleElements.h"
#include "QuadElements.h"
#include "TetrahedronElements.h"
#include "HexahedronElements.h"

class Data {
public:
	Data() {}
	Data(std::shared_ptr <Parser> p);
	virtual ~Data() = default;

	const std::shared_ptr<Parser>& getParser() const { return parser; }
	const std::shared_ptr<Element> getElem(int i) const { return elements[i]; }

	Eigen::SparseMatrix <double> K;
	Eigen::SparseVector <double> R;
	Eigen::VectorX <double> U; // displacement

	void solve();

private:
	std::vector <shared_ptr<Element>> elements;
	std::vector <Node> nodes;

	std::shared_ptr<Parser> parser;
	int dim;

	void create_nodes();
	void create_elements();
	void create_constrains();
	void create_load();

	void fillGlobalK();
	void fillGlobalR();

	void fillConstrains();
	void addToGlobalK(int first_index, int second_index, double value);
	void addToGlobalR(int index, double value);
};


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
