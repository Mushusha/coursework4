#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Core"

#include "Element.h"
#include "MathMV.h"

class quadElement : public Element {
public:
	quadElement(int id, int type, std::vector <int> nodes)
		: Element(id, type, nodes) {}
	virtual ~quadElement() = default;

	Eigen::MatrixXd localK() override;
	std::vector <double> localF() override;
	Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) override;
	std::vector<double> FF(double ksi, double eta, double zeta = 0) override;

	bool pointInElem(std::vector<double> point) override;

	Eigen::MatrixXd localC() override;
	std::vector <double> localR(std::vector<double> value) override;
	Eigen::MatrixXd localM() override;

	std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;

protected:
	quadElement() : Element() {}

	Eigen::MatrixXd  gradFF(double ksi, double eta, double zeta = 0) override;
	Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) override;

	void set_pressure(int edge, double value) override;
private:
	double len_edge(int edge);
	std::pair<int, int> edge_to_node(int edge);
};


class infQuadElement : public  quadElement {
public:
	infQuadElement(int id, int type, std::vector<int> nodes)
		: quadElement(id, type, nodes) {
	}
	virtual ~infQuadElement() = default;

	std::vector<double> FF(double ksi, double eta, double zeta = 0);

private:
	infQuadElement() : quadElement() {}
};
