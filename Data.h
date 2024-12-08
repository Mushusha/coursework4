#pragma once

#include "Parser.h"
#include "MathMV.h"
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

	const std::shared_ptr<Parser>& get_parser() const { return parser; }
	const std::shared_ptr<Element> get_elem(int i) const { return elements[i]; }

	Eigen::SparseMatrix <double> K;
	Eigen::SparseVector <double> F;
	Eigen::VectorX <double> U; // displacement

	Eigen::SparseMatrix <double> C; // agreed resultants
	Eigen::SparseVector <double> R;


	void fillGlobalK();
	void fillGlobalF();
	void fillconstraints();

	void fillGlobalC();
	void fillGlobalR(std::string field);

	void solve();
	void fillFields();
	void smoothing(std::string field);

private:
	std::vector <shared_ptr<Element>> elements;
	std::vector <Node> nodes;

	std::shared_ptr<Parser> parser;
	int dim;

	void create_nodes();
	void create_elements();
	void create_constants();
	void create_constraints();
	void create_load();

	void addToGlobalK(int first_index, int second_index, double value);
	void addToGlobalF(int index, double value);

	void displacementToElements();
	void displacementToNodes();
	void calcStrain();
	void calcStress();

	void zeroDiagonalCheck();

	void outputValues(std::string field);
};