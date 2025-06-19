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
#include "Node/Node.h"
#include "log1.h"

class Element {
public:
	Element(int id, ElemType type, std::vector <int> nodes);
	Element(const Element& other);
	Element& operator=(const Element& other);
	Element(Element&& other) noexcept;
	Element& operator=(Element&& other) noexcept; 
	virtual ~Element() = default;

	virtual std::vector <std::complex<double>> FF(double ksi, double eta, double zeta = 0) = 0;
	virtual Eigen::MatrixXcd B(double ksi = 0, double eta = 0, double zeta = 0) = 0;

	virtual Eigen::MatrixXcd localK() = 0; // <NODES * DIM, NODES * DIM>
	virtual std::vector <double> localF(double mult = 1) = 0;

	virtual Eigen::MatrixXd localC() = 0; // agreed resultants
	virtual std::vector <double> localR(std::vector<double> value) = 0;
	virtual Eigen::MatrixXcd localM() = 0;

	std::complex<double> Jac(double ksi, double eta, double zeta = 0);
	virtual double gaussPoint(LocVar var, int i) = 0;
	virtual double weight(LocVar var, int i) = 0;
	virtual double Volume() = 0;

	virtual std::vector<int> edge_to_node(int edge) = 0;

	Eigen::MatrixXd D;

	void set_coords(std::vector <double> x, std::vector <double> y, std::vector <double> z);
	void set_load(int type, int edge, std::array<double, 6> value); // nodes
	void set_constants(double E, double nu, double rho);

	int get_node(int i) const { return nodes[i]; }
	std::vector<int> get_node() const { return nodes; }
	int nodes_count() const { return nodes.size(); }
	double get_coord(int loc_node, int comp) const;
	int get_id() const { return id; }
	ElemType get_type() const { return type; }

	double get_E() const { return Young; }
	double get_nu() const { return Poisson; }
	double get_rho() const { return density; }
	
	std::vector<std::map <ResType, Eigen::VectorXcd>> results; // [i] - node; [j] - field; [k] - comp // tensor
	Eigen::VectorXd displacements;

	// for interpolation
	virtual bool pointInElem(std::vector<double> point) = 0;
	virtual std::vector<double> coordFF(double x0, double y0, double z0 = 0) = 0;

protected:
	Element() = default;

	ElemType type;
	int id;
	std::vector <double> x, y, z;
	std::vector <int> nodes; // find in Node

	double Young;
	double Poisson;
	double density;
	
	// refactoring: may be <edge, vector>
	std::map <std::pair<int, int>, double> load; // pair <edge, comp>, value
	
	virtual Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta = 0) = 0;
	virtual Eigen::MatrixXcd J(double ksi, double eta, double zeta = 0) = 0; // <DIM, DIM>
	
	virtual void set_pressure(int edge, double value) = 0;
};