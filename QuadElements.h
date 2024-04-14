#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "global.h"

//dim = 2;
//numNodes = 4;

class quadElement : public Element {
public:
	quadElement() {};
	quadElement(std::vector <int> _nodes, std::vector <float> _x, std::vector <float> _y);
	virtual ~quadElement() = default;
};

class quadInfElement : public quadElement {
public:
	quadInfElement() {};
	quadInfElement(std::vector <int> _nodes, std::vector <float> _x, std::vector <float> _y);
	virtual ~quadInfElement() = default;
};

float quadN1(float ksi, float eta);
float quadN2(float ksi, float eta);
float quadN3(float ksi, float eta);
float quadN4(float ksi, float eta);

float quadInfN(float ksi, float eta, std::vector <float> a);
float quadMapN(float ksi, float eta, std::vector <float> a);
//float shapeN(float ksi, float eta)

Eigen::Matrix2f quadMatrixJ(float ksi, float eta, std::vector <float> x, std::vector <float> y, float (*N)(float, float, std::vector <float>));
float quadDksiN(float ksi, float eta, float (*N)(float, float));
float quadDetaN(float ksi, float eta, float (*N)(float, float));

typedef float (*arrFunc) (float ksi, float eta);

Eigen::MatrixXf  quadMatrixB(float ksi, float eta, std::vector <float> x, std::vector <float> y, float (*funcN)(float, float, std::vector <float>));
Eigen::MatrixXf quadLocK(std::vector <float> x, std::vector <float> y, float (*N)(float, float, std::vector <float>));