#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "global.h"

//dim = 2;
//numNodes = 3;

class triangleElement : public Element {
public:
	triangleElement() {};
	triangleElement(std::vector <int> _nodes, std::vector <float> _x, std::vector <float> _y);
	virtual ~triangleElement() = default;
};

Eigen::Matrix3f triangleMatrixC(std::vector <float> x, std::vector <float> y);
Eigen::MatrixXf triangleMatrixB(std::vector <float> x, std::vector <float> y);
Eigen::MatrixXf triangleLocK(std::vector <float> x, std::vector <float> y);