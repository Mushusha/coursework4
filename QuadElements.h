#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class quadElement : public Element {
public:
	quadElement() : Element() {}
	quadElement(int id, int type, std::vector <int> nodes)
		: Element(id, type, nodes) {}
	virtual ~quadElement() = default;

	Eigen::MatrixXd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(double value) override;

protected:
	std::vector<double> FF(double ksi, double eta, double zeta = 0) override;
	std::vector <std::vector <double>> gradFF(double ksi, double eta, double zeta = 0) override;
	Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) override;

	void set_pressure(int edge, double value);

};

//class quadInfElement : public quadElement {
//public:
//	quadInfElement() {};
//	quadInfElement(std::vector <int> _nodes, std::vector <double> _x, std::vector <double> _y);
//	virtual ~quadInfElement() = default;
//};