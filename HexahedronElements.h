#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class hexElement : public Element {
public:
	hexElement() {};
	hexElement(int id, int type, std::vector <int> nodes)
		: Element (id, type, nodes) {};
	virtual ~hexElement() = default;

	Eigen::MatrixXd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXd B(double ksi, double eta, double zeta) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(double value) override;

protected:
	std::vector<double> FF(double ksi, double eta, double zeta) override;
	std::vector<std::vector<double>> gradFF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXd J(double ksi, double eta, double zeta) override;

	void set_pressure(int edge, double value);

};

//class hexaInfElement : public hexElement {
//public:
//	hexaInfElement() {};
//	virtual ~hexaInfElement() = default;
//	/*
//	dim = 3;
//	numNodes = 8;
//	quadLocK(hexaInfN);
//	*/
//};

