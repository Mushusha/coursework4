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
#include "Eigen/Core"

#include "Parser.h"
#include "Enums.h"
#include "Node.h"


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


	virtual std::vector <double> FF(double ksi, double eta, double zeta = 0) = 0;
	virtual std::vector <std::vector <double>> gradFF(double ksi, double eta, double zeta = 0) = 0;
	virtual Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) = 0;
	virtual Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) = 0; // <DIM, DIM>
	virtual Eigen::MatrixXd locK() = 0; // <NODES * DIM, NODES * DIM>
	virtual std::vector <double> locR() = 0;


	void set_coords(std::vector <double> x, std::vector <double> y, std::vector <double> z);
	void set_load(int edge, int comp, double value);

	int get_nodes(int i) { return nodes[i]; }
	int nodes_count() { return nodes.size(); }
protected:
	int type;
	int id;
	std::vector <double> x, y, z;
	std::vector <int> nodes; // find in Node

	std::map <std::pair<int, int>, double> load; // pair <loc edge, comp>, value


	Eigen::MatrixXd twoMatrixD(); //plane_stress, plane_strain ??  
	Eigen::MatrixXd threeMatrixD();
};

//Eigen::SparseMatrix <double> globalK(std::shared_ptr <Parser> p, int numNodes, int dim);
//Eigen::SparseVector <double> globalLoad(std::shared_ptr <Parser> p, int numNodes, int dim);
//int indexConstrain(std::shared_ptr <Parser> p, int count, int apply, int dim);
