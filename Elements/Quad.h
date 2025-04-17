#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

#include "Eigen/Core"

#include "Element.h"
#include "MathMV.h"

class Quad : public Element {
public:
	Quad() : Element() {}
	Quad(int id, ElemType type, std::vector <int> nodes)
		: Element(id, type, nodes) {}
	Quad(const Quad& other)
		: Element(other) {}
	Quad& operator=(const Quad& other) {
		if (this != &other)
			Element::operator=(other);
		return *this;
	}
	Quad(Quad&& other) noexcept
		: Element(std::move(other)) {}
	Quad& operator=(Quad&& other) noexcept {
		if (this != &other)
			Element::operator=(std::move(other));
		return *this;
	}
	virtual ~Quad() = default;

	Eigen::MatrixXcd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXcd B(double ksi = 0, double eta = 0, double zeta = 0) override;
	std::vector<std::complex<double>> FF(double ksi, double eta, double zeta = 0) override;

	bool pointInElem(std::vector<double> point) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(std::vector<double> value) override;
	Eigen::MatrixXcd localM() override;

	std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;
	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;
	double Volume() final;

protected:
	Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta = 0) override;
	Eigen::MatrixXcd J(double ksi, double eta, double zeta = 0) override;

	void set_pressure(int edge, double value) override;
private:
	double len_edge(int edge);
	std::pair<int, int> edge_to_node(int edge);
};


class infQuad : public  Quad {
public:
	infQuad() : Quad() {}
	infQuad(int id, ElemType type, std::vector<int> nodes)
		: Quad(id, type, nodes) {
	}
	infQuad(const infQuad& other)
		: Quad(other) {}
	infQuad& operator=(const infQuad& other) {
		if (this != &other)
			Element::operator=(other);
		return *this;
	}
	infQuad(infQuad&& other) noexcept
		: Quad(std::move(other)) {}
	infQuad& operator=(infQuad&& other) noexcept {
		if (this != &other)
			Element::operator=(std::move(other));
		return *this;
	}
	virtual ~infQuad() = default;

	bool is_dyn;
	double omega;

	std::vector<std::complex<double>> FF(double ksi, double eta, double zeta = 0) final;
	double gaussPoint(LocVar var, int i) final;
	double weight(LocVar var, int i) final;

private:
	double A;
	double k;
};