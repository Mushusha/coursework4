#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class Wedge : public Element {
public:
	Wedge() : Element() {};
	Wedge(int id, ElemType type, std::vector <int> nodes)
		: Element(id, type, nodes) {
	};
	Wedge(const Wedge& other)
		: Element(other) {
	}
	Wedge& operator=(const Wedge& other) {
		if (this != &other)
			Element::operator=(other);
		return *this;
	}
	Wedge(Wedge&& other) noexcept
		: Element(std::move(other)) {
	}
	Wedge& operator=(Wedge&& other) noexcept {
		if (this != &other) {
			Element::operator=(std::move(other));
		}
		return *this;
	}
	virtual ~Wedge() = default;

	std::vector<std::complex<double>> FF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXcd B(double ksi, double eta, double zeta) override;

	Eigen::MatrixXcd localK() override;
	std::vector <double> localF(double mult = 1) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(std::vector<double> value) override;
	Eigen::MatrixXcd localM() override;

	std::vector<int> edge_to_node(int edge) final;

	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;
	double Volume() final;

	bool pointInElem(std::vector<double> point) override;
	std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;

protected:
	Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXcd J(double ksi, double eta, double zeta) override;

	void set_pressure(int edge, double value);
	std::array<double, 3> normal(int edge);
	double area_edge(int edge);
};

