#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

#include "Eigen/Core"

#include "Quad.h"


/**
 * Infinite Quad Element (2D)
 * 
 * Standard mapped infinite element for 2D analysis.
 * Uses decay functions in the radial (ksi) direction.
 * 
 * Node numbering (4 nodes):
 *   Boundary at ksi = -1 (interface with finite element): nodes 0, 3
 *   At infinity ksi = 1: nodes 1, 2
 * 
 *     3-----------2     (eta = 1)
 *     |           |
 *     |           |
 *     |           |
 *     0-----------1     (eta = -1)
 *   ksi=-1      ksi=1
 *   (boundary)  (infinity)
 */
class InfQuad : public Quad {
public:
	InfQuad() : Quad(), is_dyn(false), omega(0.0) {}
	InfQuad(int id, ElemType type, std::vector<int> nodes)
		: Quad(id, type, nodes), is_dyn(false), omega(0.0) {}
	InfQuad(const InfQuad& other)
		: Quad(other), is_dyn(other.is_dyn), omega(other.omega) {}
	InfQuad& operator=(const InfQuad& other) {
		if (this != &other) {
			Quad::operator=(other);
			is_dyn = other.is_dyn;
			omega = other.omega;
		}
		return *this;
	}
	InfQuad(InfQuad&& other) noexcept
		: Quad(std::move(other)), is_dyn(other.is_dyn), omega(other.omega) {}
	InfQuad& operator=(InfQuad&& other) noexcept {
		if (this != &other) {
			Quad::operator=(std::move(other));
			is_dyn = other.is_dyn;
			omega = other.omega;
		}
		return *this;
	}
	virtual ~InfQuad() = default;

	bool is_dyn = false;
	double omega = 0.0;
	
	double pole_x = 0.0, pole_y = 0.0, pole_z = 0.0;
	void set_pole(double px, double py, double pz = 0.0) {
		pole_x = px; pole_y = py; pole_z = pz;
	}

	std::vector<std::complex<double>> FF(double ksi, double eta, double zeta = 0) override;
	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;

protected:
	Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta = 0) override;

private:
	double compute_A() const;
};