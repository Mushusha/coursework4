#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class triElement : public Element {
public:
	triElement() {};
	triElement(int id, int type, std::vector <int> nodes)
		: Element (id, type, nodes) {};
	virtual ~triElement() = default;
	
	Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) override;
	Eigen::MatrixXd locK() override;

	std::vector <double> FF(double ksi, double eta, double zeta = 0) override;
	std::vector <std::vector <double>> gradFF(double ksi, double eta, double zeta = 0) override;
	Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) override;

private:
	Eigen::MatrixXd C();
};

