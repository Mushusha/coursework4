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
#include "log.h"

class Element {
public:
	Element() {};
	Element(int id, int type, std::vector <int> nodes) {
		this->id = id;
		this->type = type;
		this->nodes.resize(nodes.size());
		for (int i = 0; i < nodes.size(); i++)
			this->nodes[i] = nodes[i];
	};
	virtual ~Element() = default;

	virtual Eigen::MatrixXd localK() = 0; // <NODES * DIM, NODES * DIM>
	virtual std::vector <double> localF() = 0;
	virtual Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) = 0;
	virtual std::vector <double> FF(double ksi, double eta, double zeta = 0) = 0;

	virtual bool pointInElem(std::vector<double> point) = 0;

	virtual Eigen::MatrixXd localC() = 0; // agreed resultants
	virtual std::vector <double> localR(std::vector<double> value) = 0;

	Eigen::MatrixXd planeStressD(); //plane_stress, plane_strain ??  
	Eigen::MatrixXd planeStrainD();
	Eigen::MatrixXd threeMatrixD();

	void set_coords(std::vector <double> x, std::vector <double> y, std::vector <double> z);
	void set_load(int type, int edge, std::array<double, 6> value); // nodes
	void set_constants(double E, double nu);

	int get_type() { return type; }
	int get_nodes(int i) { return nodes[i]; }
	int nodes_count() { return nodes.size(); }
	double get_coord(int loc_node, int comp);

	std::vector<std::vector <Eigen::VectorXd>> results; // [i] - node; [j] - field; [k] - comp
	Eigen::VectorXd displacements;
protected:
	int type;
	int id;
	std::vector <double> x, y, z;
	std::vector <int> nodes; // find in Node

	double Young;
	double Poisson;

	// refactoring: may be <edge, vector>
	std::map <std::pair<int, int>, double> load; // pair <edge, comp>, value

	virtual Eigen::MatrixXd gradFF(double ksi, double eta, double zeta = 0) = 0;
	virtual Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) = 0; // <DIM, DIM>
	
	virtual void set_pressure(int edge, double value) = 0;
};