#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class tetraElement : public Element {
public:
	tetraElement() : Element() {};
	tetraElement(int id, int type, std::vector <int> nodes)
		: Element (id, type, nodes) {};
	tetraElement(const tetraElement& other)
		: Element(other) {}
	tetraElement& operator=(const tetraElement& other) {
		if (this != &other)
			Element::operator=(other);
		return *this;
	}
	tetraElement(tetraElement&& other) noexcept
		: Element(std::move(other)) {}
	tetraElement& operator=(tetraElement&& other) noexcept {
		if (this != &other)
			Element::operator=(std::move(other));
		return *this;
	}
	virtual ~tetraElement() = default;

	Eigen::MatrixXd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) override;
	std::vector <double> FF(double ksi, double eta, double zeta) override;

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