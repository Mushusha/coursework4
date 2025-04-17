#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class Tetra : public Element {
public:
	Tetra() : Element() {};
	Tetra(int id, ElemType type, std::vector <int> nodes)
		: Element (id, type, nodes) {};
	Tetra(const Tetra& other)
		: Element(other) {}
	Tetra& operator=(const Tetra& other) {
		if (this != &other)
			Element::operator=(other);
		return *this;
	}
	Tetra(Tetra&& other) noexcept
		: Element(std::move(other)) {}
	Tetra& operator=(Tetra&& other) noexcept {
		if (this != &other)
			Element::operator=(std::move(other));
		return *this;
	}
	virtual ~Tetra() = default;

	Eigen::MatrixXcd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXcd B(double ksi = 0, double eta = 0, double zeta = 0) override;
	std::vector <std::complex<double>> FF(double ksi, double eta, double zeta) override;

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

private:
	Eigen::MatrixXd C();
};

//class tetraInfElement : public Tetra {
//public:
//	tetraInfElement() {};
//	virtual ~tetraInfElement() = default;
//	/*
//	dim = 3;
//	numNodes = 4;
//	tetraLocK();
//	*/
//};