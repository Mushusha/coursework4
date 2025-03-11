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
	Data(const Data& other)
		: dim(other.dim),
		line_start(other.line_start),
		line_end(other.line_end),
		points_count(other.points_count),
		out_stress(other.out_stress),
		elements(other.elements),
		nodes(other.nodes),
		parser(other.parser) {}
	Data& operator=(const Data& other) {
		if (this != &other) {
			dim = other.dim;
			line_start = other.line_start;
			line_end = other.line_end;
			points_count = other.points_count;
			out_stress = other.out_stress;
			elements = other.elements;
			nodes = other.nodes;
			parser = other.parser;
		}
		return *this;
	}
	Data(Data&& other) noexcept
		: dim(std::exchange(other.dim, 0)),
		line_start(std::move(other.line_start)),
		line_end(std::move(other.line_end)),
		points_count(std::exchange(other.points_count, 0)),
		out_stress(std::move(other.out_stress)),
		elements(std::move(other.elements)),
		nodes(std::move(other.nodes)),
		parser(std::move(other.parser)) {}
	Data& operator=(Data&& other) noexcept {
		if (this != &other) {
			dim = std::exchange(other.dim, 0);
			line_start = std::move(other.line_start);
			line_end = std::move(other.line_end);
			points_count = std::exchange(other.points_count, 0);
			out_stress = std::move(other.out_stress);
			elements = std::move(other.elements);
			nodes = std::move(other.nodes);
			parser = std::move(other.parser);
		}
		return *this;
	}
	virtual ~Data() = default;

	const std::shared_ptr<Parser>& get_parser() const { return parser; }
	const std::shared_ptr<Element> get_elem(int i) const { return elements[i]; }
	void set_output_param(std::vector<double> start, std::vector<double> end, int count);

	Eigen::SparseMatrix <double> K;
	Eigen::SparseVector <double> F;
	Eigen::SparseMatrix <double> M;
	Eigen::VectorX <double> U; // displacement

	int dim;

	std::vector<double> line_start;
	std::vector<double> line_end;
	int points_count;

	std::vector<std::vector<double>> out_stress;

	void fillGlobalK();
	void fillGlobalF();
	void fillGlobalM();
	void fillConstraints();

	void fillGlobalC();
	void fillGlobalR(int type, int comp);

	void solve();

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

	void addToGlobalK(int first_index, int second_index, double value);
	void addToGlobalF(int index, double value);

	void displacementToElements();
	void displacementToNodes();
	void calcStrain();
	void calcStress();

	void fillFields();
	void smoothing();

	void zeroDiagonalCheck();

	void outputValues(int type, int comp);
	void printMeshStress();

	void interpolation(std::vector<std::vector<double>>& points, std::vector<double>& values, int type, int comp);

};