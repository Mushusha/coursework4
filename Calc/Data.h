#pragma once

#include "Parser.h"
#include "log1.h"
#include "Enums.h"
#include "Node/Node.h"
#include "Element.h"
#include "Tri.h"
#include "Quad.h"
#include "Tetra.h"
#include "Hex.h"


class Data {
public:
	Data(std::shared_ptr<const Parser> p);
	Data(const Data& other);
	Data& operator=(const Data& other);
	Data(Data&& other) noexcept;
	Data& operator=(Data&& other) noexcept;
	virtual ~Data() = default;

	const std::shared_ptr<Element> get_elem(int i) const { return elements[i]; }
	const std::shared_ptr<Node> get_node(int i) const { return nodes[i]; }
	const int nodes_count() const { return nodes.size(); }
	const int elements_count() const { return elements.size(); }

	std::vector<shared_ptr<Element>> get_elements() { return elements; }
	std::vector<shared_ptr<Node>> get_nodes() { return nodes; }

	int dim;
	std::string analisys_type;

	double damping;
	double max_time;
	int max_iter;
	int num_inf_elems;
	double omega;
	double Amp;
	  
private:
	Data() = default;

	std::vector<shared_ptr<Element>> elements;
	std::vector<shared_ptr<Node>> nodes;

	void create_nodes(std::shared_ptr <const Parser> parser);
	void create_elements(std::shared_ptr <const Parser> parser);
	void create_infelements(std::shared_ptr <const Parser> parser);
	void create_constants(std::shared_ptr <const Parser> parser);
	void create_constraints(std::shared_ptr <const Parser> parser);
	void create_load(std::shared_ptr <const Parser> parser);
	void create_D(std::shared_ptr <const Parser> parser);
}; 