#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class Hex : public Element {
public:
	Hex() : Element() {};
	Hex(int id, int type, std::vector <int> nodes)
		: Element (id, type, nodes) {};
	Hex(const Hex& other)
		: Element(other) {}
	Hex& operator=(const Hex& other) {
		if (this != &other)
			Element::operator=(other);
		return *this;
	}
	Hex(Hex&& other) noexcept
		: Element(std::move(other)) {}
	Hex& operator=(Hex&& other) noexcept {
		if (this != &other) {
			Element::operator=(std::move(other));
		}
		return *this;
	}
	virtual ~Hex() = default;

	Eigen::MatrixXd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXd B(double ksi, double eta, double zeta) override;
	std::vector<double> FF(double ksi, double eta, double zeta) override;

	bool pointInElem(std::vector<double> point) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(std::vector<double> value) override;
	Eigen::MatrixXd localM() override;

	std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;

	double Volume() final;

protected:
	Eigen::MatrixXd gradFF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXd J(double ksi, double eta, double zeta) override;

	void set_pressure(int edge, double value);
};

//class hexaInfElement : public Hex {
//public:
//	hexaInfElement() {};
//	virtual ~hexaInfElement() = default;
//	/*
//	dim = 3;
//	numNodes = 8;
//	quadLocK(hexaInfN);
//	*/
//};

