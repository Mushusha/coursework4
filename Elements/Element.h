#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <cstring>
#include <cmath>
#include <memory>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/Core"
#include "unsupported/Eigen/IterativeSolvers"

#include "Parser.h"
#include "Enums.h"
#include "Node.h"
#include "log1.h"

class Element {
public:
	Element(int id, int type, std::vector <int> nodes) :
		id(id), type(type), nodes(nodes) {
	};
	Element(const Element& other)
		: id(other.id),
		type(other.type),
		nodes(other.nodes),
		x(other.x),
		y(other.y),
		z(other.z),
		Young(other.Young),
		Poisson(other.Poisson),
		density(other.density),
		load(other.load),
		results(other.results),
		displacements(other.displacements),
		D(other.D) {}
	Element& operator=(const Element& other) {
		if (this != &other) {
			id = other.id;
			type = other.type;
			nodes = other.nodes;
			x = other.x;
			y = other.y;
			z = other.z;
			Young = other.Young;
			Poisson = other.Poisson;
			density = other.density;
			load = other.load;
			results = other.results;
			displacements = other.displacements;
			D = other.D;
		}
		return *this;
	}
	Element(Element&& other) noexcept
		: id(std::exchange(other.id, 0)),
		type(std::exchange(other.type, 0)),
		nodes(std::move(other.nodes)),
		x(std::move(other.x)),
		y(std::move(other.y)),
		z(std::move(other.z)),
		Young(std::exchange(other.Young, 0.0)),
		Poisson(std::exchange(other.Poisson, 0.0)),
		density(std::exchange(other.density, 0.0)),
		load(std::move(other.load)),
		results(std::move(other.results)),
		displacements(std::move(other.displacements)),
		D(std::move(other.D)) {
	}
	Element& operator=(Element&& other) noexcept {
		if (this != &other) {
			id = std::exchange(other.id, 0);
			type = std::exchange(other.type, 0);
			nodes = std::move(other.nodes);
			x = std::move(other.x);
			y = std::move(other.y);
			z = std::move(other.z);
			Young = std::exchange(other.Young, 0.0);
			Poisson = std::exchange(other.Poisson, 0.0);
			density = std::exchange(other.density, 0.0);
			load = std::move(other.load);
			results = std::move(other.results);
			displacements = std::move(other.displacements);
			D = std::move(other.D);
		}
		return *this;
	}
	virtual ~Element() = default;

	virtual Eigen::MatrixXd localK() = 0; // <NODES * DIM, NODES * DIM>
	virtual std::vector <double> localF() = 0;
	virtual Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) = 0;
	virtual std::vector <double> FF(double ksi, double eta, double zeta = 0) = 0;

	virtual bool pointInElem(std::vector<double> point) = 0;

	virtual Eigen::MatrixXd localC() = 0; // agreed resultants
	virtual std::vector <double> localR(std::vector<double> value) = 0;
	virtual Eigen::MatrixXd localM() = 0;

	virtual std::vector<double> coordFF(double x0, double y0, double z0 = 0) = 0;

	Eigen::MatrixXd D;

	void set_coords(std::vector <double> x, std::vector <double> y, std::vector <double> z);
	void set_load(int type, int edge, std::array<double, 6> value); // nodes
	void set_constants(double E, double nu, double rho);

	int get_type() { return type; }
	int get_nodes(int i) { return nodes[i]; }
	int nodes_count() { return nodes.size(); }
	double get_coord(int loc_node, int comp);

	double get_E() { return Young; }
	double get_nu() { return Poisson; }
	
	double Jac(double ksi, double eta, double zeta = 0);
	
	std::vector<std::vector <Eigen::VectorXd>> results; // [i] - node; [j] - field; [k] - comp
	Eigen::VectorXd displacements;

protected:
	Element() = default;

	int type;
	int id;
	std::vector <double> x, y, z;
	std::vector <int> nodes; // find in Node

	double Young;
	double Poisson;
	double density;
	
	// refactoring: may be <edge, vector>
	std::map <std::pair<int, int>, double> load; // pair <edge, comp>, value
	
	virtual Eigen::MatrixXd gradFF(double ksi, double eta, double zeta = 0) = 0;
	virtual Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) = 0; // <DIM, DIM>
	
	virtual void set_pressure(int edge, double value) = 0;
};