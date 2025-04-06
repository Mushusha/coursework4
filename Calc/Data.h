#pragma once

#include "Parser.h"
#include "MathMV.h"
#include "log1.h"
#include "Enums.h"
#include "Node.h"
#include "Element.h"
#include "Tri.h"
#include "Quad.h"
#include "Tetra.h"
#include "Hex.h"


class Data {
public:
	Data(std::shared_ptr<Parser> p);
	Data(const Data& other);
	Data& operator=(const Data& other);
	Data(Data&& other) noexcept;
	Data& operator=(Data&& other) noexcept;
	virtual ~Data() = default;

	const std::shared_ptr<Parser>& get_parser() const { return parser; }
	const std::shared_ptr<Element> get_elem(int i) const { return elements[i]; }
	const std::shared_ptr<Node> get_node(int i) const { return nodes[i]; }
	const int nodes_count() const { return nodes.size(); }
	const int elements_count() const { return elements.size(); }

	int dim;
	std::string analisys_type;

	std::vector<std::vector<double>> out_stress;

	double damping;
	double max_time;
	int max_iter;
	  
private:
	Data() = default;

	std::vector<shared_ptr<Element>> elements;
	std::vector<shared_ptr<Node>> nodes;

	std::shared_ptr<Parser> parser;

	void create_nodes();
	void create_elements();
	void create_infelements();
	void create_constants();
	void create_constraints();
	void create_load();
	void create_D();
}; 