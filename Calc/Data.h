#pragma once

#include "Parser.h"
#include "MathMV.h"
#include "log1.h"
#include "Enums.h"
#include "Node.h"
#include "Element.h"
#include "TriangleElements.h"
#include "QuadElements.h"
#include "TetrahedronElements.h"
#include "HexahedronElements.h"


class Data {
public:
	Data(std::shared_ptr <Parser> p);
	Data(const Data& other);
	Data& operator=(const Data& other);
	Data(Data&& other) noexcept;
	Data& operator=(Data&& other) noexcept;
	virtual ~Data() = default;

	const std::shared_ptr<Parser>& get_parser() const { return parser; }
	const std::shared_ptr<Element> get_elem(int i) const { return elements[i]; }
	Node get_node(int i) const { return nodes[i]; }
	const int nodes_count() const { return nodes.size(); }
	const int elements_count() const { return elements.size(); }
	void set_output_param(std::vector<double> start, std::vector<double> end, int count);

	int dim;
	std::string analisys_type;

	std::vector<double> line_start;
	std::vector<double> line_end;
	int points_count;

	std::vector<std::vector<double>> out_stress;

	//void fillGlobalC();
	//void fillGlobalR(int type, int comp);

private:
	Data() = default;

	std::vector <shared_ptr<Element>> elements;
	std::vector <Node> nodes;

	std::shared_ptr<Parser> parser;

	Eigen::SparseMatrix <double> C; // agreed resultants
	Eigen::SparseVector <double> R;

	void create_nodes();
	void create_elements();
	void create_infelements();
	void create_constants();
	void create_constraints();
	void create_load();
	void create_D();

	//void smoothing();

	//void outputValues(int type, int comp);
	//void printMeshStress();

	//void interpolation(std::vector<std::vector<double>>& points, std::vector<double>& values, int type, int comp);

};