#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "global.h"

class tetraElement : public Element {
public:
	tetraElement() {};
	virtual ~tetraElement() = default;
	/*
	dim = 3;
	numNodes = 4;
	tetraLocK();
	*/
};

class tetraInfElement : public tetraElement {
public:
	tetraInfElement() {};
	virtual ~tetraInfElement() = default;
	/*
	dim = 3;
	numNodes = 4;
	tetraLocK();
	*/
};

Eigen::Matrix4f tetraMatrixC(std::vector <float> x, std::vector <float> y, std:: vector <float> z);
Eigen::MatrixXf tetraMatrixB(std::vector <float> x, std::vector <float> y, std::vector <float> z);
Eigen::MatrixXf tetraLocK(std::vector <float> x, std::vector <float> y, std::vector <float> z);