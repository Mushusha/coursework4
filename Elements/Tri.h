#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"


class Tri : public Element {
public:
	Tri() : Element() {};
	Tri(int id, ElemType type, std::vector <int> nodes)
		: Element (id, type, nodes) {};
	Tri(const Tri& other)
		: Element(other) {}
	Tri& operator=(const Tri& other) {
		if (this != &other)
			Element::operator=(other);
		return *this;
	}
	Tri(Tri&& other) noexcept
		: Element(std::move(other)) {}
	Tri& operator=(Tri&& other) noexcept {
		if (this != &other)
			Element::operator=(std::move(other));
		return *this;
	}
	virtual ~Tri() = default;
	
	Eigen::MatrixXd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) override;
	std::vector <double> FF(double ksi, double eta, double zeta = 0) override;

	bool pointInElem(std::vector<double> point) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(std::vector<double> value) override;
	Eigen::MatrixXd localM() override;

	std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;
	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;
	double Volume() final;

protected:
	Eigen::MatrixXd gradFF(double ksi, double eta, double zeta = 0) override;
	Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) override;

	void set_pressure(int edge, double value);

private:
	Eigen::MatrixXd C();
	double len_edge(int edge);
	std::pair<int, int> edge_to_node(int edge);
};

