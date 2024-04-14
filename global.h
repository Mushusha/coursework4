#pragma once
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <memory>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"

#include "Parser.h"

class Element {
public:
	Element() {};
	virtual ~Element() = default;
	Eigen::MatrixXf locK;
protected:
	std::vector <int> nodes;
	std::vector <float> x;
	std::vector <float> y;
	std::vector <float> z;
};

Eigen::Matrix3f twoMatrixD();
Eigen::MatrixXf threeMatrixD();

Eigen::SparseMatrix <float> globalK(std::shared_ptr <Parser> p, int numNodes, int dim);

std::vector <std::shared_ptr <Element>> createElements(std::shared_ptr <Parser> p, int dim);