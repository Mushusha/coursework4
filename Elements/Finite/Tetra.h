#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"
#include "Data.h"


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

	std::vector <std::complex<double>> FF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXcd B(double ksi = 0, double eta = 0, double zeta = 0) override;

	Eigen::MatrixXcd localK() override;
	std::vector <double> localF(double mult = 1) override;

	std::vector<int> edge_to_node(int edge) final;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(std::vector<double> value) override;
	Eigen::MatrixXcd localM() override;
	Eigen::MatrixXd localDamping() override;

	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;
	double Volume() final;

	bool pointInElem(std::vector<double> point) override;
	std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;

protected:
	Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXcd J(double ksi, double eta, double zeta) override;

	void set_pressure(int edge, double value);

private:
	Eigen::MatrixXd C();

	std::array<double, 3> normal(int edge);
	double area_edge(int edge);
};