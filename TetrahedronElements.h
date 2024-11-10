#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class tetraElement : public Element {
public:
	tetraElement() {};
	tetraElement(int id, int type, std::vector <int> nodes)
		: Element (id, type, nodes) {};

	virtual ~tetraElement() = default;

	Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) override;
	Eigen::MatrixXd locK() override;
	std::vector <double> locR() override;

	std::vector <double> FF(double ksi, double eta, double zeta) override;
	std::vector <std::vector <double>> gradFF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXd J(double ksi, double eta, double zeta) override;

private:
	Eigen::MatrixXd C();
};

//class tetraInfElement : public tetraElement {
//public:
//	tetraInfElement() {};
//	virtual ~tetraInfElement() = default;
//	/*
//	dim = 3;
//	numNodes = 4;
//	tetraLocK();
//	*/
//};