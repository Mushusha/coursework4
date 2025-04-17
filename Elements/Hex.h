#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class Hex : public Element {
public:
	Hex() : Element() {};
	Hex(int id, ElemType type, std::vector <int> nodes)
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

	Eigen::MatrixXcd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXcd B(double ksi, double eta, double zeta) override;
	std::vector<std::complex<double>> FF(double ksi, double eta, double zeta) override;

	bool pointInElem(std::vector<double> point) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(std::vector<double> value) override;
	Eigen::MatrixXcd localM() override;

	std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;
	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;
	double Volume() final;

protected:
	Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXcd J(double ksi, double eta, double zeta) override;

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

