#pragma once

#include <cmath>
#include <algorithm>
#include <map>
#include <execution>
#include <numeric>

#include "Parser.h"
#include "log1.h"
#include "Enums.h"
#include "MathMV.h"
#include "Node/Node.h"
#include "Element.h"
#include "Tri.h"
#include "Quad.h"
#include "Tetra.h"
#include "Hex.h"
#include "Wedge.h"
#include "Pyr.h"
#include "InfQuad.h"
#include "InfHex.h"
#include "QuadSEM.h"
#include "HexSEM.h"
#include "InfQuadSEM.h"
#include "InfHexSEM.h"


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

	std::string filename;

	double max_time;
	int max_iter;
	int iter_res_output;

	int num_inf_elems;

	double omega;
	double Amp;
	
	static bool is_dynamic;
	static double reference_element_size;
	static double damping;
	  
private:
	Data() = default;

	std::vector<shared_ptr<Element>> elements;
	std::vector<shared_ptr<Node>> nodes;
	std::map<int, int> node_id_old_to_new;

	void create_nodes(std::shared_ptr <const Parser> parser);
	void create_elements(std::shared_ptr <const Parser> parser);
	void create_infelements(std::shared_ptr <const Parser> parser);
	void create_constants(std::shared_ptr <const Parser> parser);
	void create_constraints(std::shared_ptr <const Parser> parser);
	void create_load(std::shared_ptr <const Parser> parser);
	void create_D(std::shared_ptr <const Parser> parser);

	void build_spectral_quad_element(std::shared_ptr<Element> elem);

	bool is_quad_node_on_constrained_edge(int i, int j, int nodes_per_side,
		const std::vector<bool>& edge_fully_constrained);
	void apply_constraints_new_quad_nodes(
		std::shared_ptr<Node> new_node, double xi, double eta,
		const std::vector<int>& original_nodes,
		const std::map<int, std::map<int, double>>& original_constraints);

	void build_spectral_hex_element(std::shared_ptr<Element> elem);
	bool is_hex_node_on_constrained_face(int i, int j, int k, 
		int nodes_per_side, const std::vector<bool>& face_fully_constrained);
	void apply_constraints_new_hex_nodes(std::shared_ptr<Node> new_node, 
		double xi, double eta, double zeta, const std::vector<int>& original_nodes, 
		const std::map<int, std::map<int, double>>& original_constraints);

	void transfer_pressure_to_infinite_element(std::shared_ptr<const Parser> parser,
		std::shared_ptr<Element> parent_elem, std::shared_ptr<Element> inf_elem,
		int parent_elem_id, int parent_side);

	void renumber_nodes();
	
	void compute_reference_element_size();
};