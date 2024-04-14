#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "global.h"

class hexaElement : public Element {
public:
	hexaElement() {};
	virtual ~hexaElement() = default;
	/*
	dim = 3;
	numNodes = 8;
	hexaLocK(hexaMapN);
	*/
};

class hexaInfElement : public hexaElement {
public:
	hexaInfElement() {};
	virtual ~hexaInfElement() = default;
	/*
	dim = 3;
	numNodes = 8;
	quadLocK(hexaInfN);
	*/
};

float hexaN1(float ksi, float eta, float zeta);
float hexaN2(float ksi, float eta, float zeta);
float hexaN3(float ksi, float eta, float zeta);
float hexaN4(float ksi, float eta, float zeta);
float hexaN5(float ksi, float eta, float zeta);
float hexaN6(float ksi, float eta, float zeta);
float hexaN7(float ksi, float eta, float zeta);
float hexaN8(float ksi, float eta, float zeta);

float hexaInfN(float ksi, float eta, float zeta, std::vector <float> a);

typedef float (*hexaArrFunc) (float ksi, float eta, float zeta);

float hexaMapN(float ksi, float eta, float zeta, std::vector <float> a);

Eigen::Matrix3f hexaMatrixJ(float ksi, float eta, float zeta, std::vector <float> x, std::vector <float> y, std::vector <float> z, float (*N)(float, float, float, std::vector <float>));

float hexaDksiN(float ksi, float eta, float zeta, float (*N)(float, float, float));
float hexaDetaN(float ksi, float eta, float zeta, float (*N)(float, float, float));
float hexaDzetaN(float ksi, float eta, float zeta, float (*N)(float, float, float));

Eigen::MatrixXf hexaMatrixB(float ksi, float eta, float zeta, std::vector <float> x, std::vector <float> y, std::vector <float> z, float (*funcN)(float, float, float, std::vector <float>));
Eigen::MatrixXf hexaLocK(std::vector <float> x, std::vector <float> y, std::vector <float> z, float (*N)(float, float, float, std::vector <float>));