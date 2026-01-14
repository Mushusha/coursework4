#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

#include "Eigen/Core"

#include "Hex.h"


class InfHex : public Hex {
public:
	InfHex() : Hex(), omega(0.0) {}
	InfHex(int id, ElemType type, std::vector<int> nodes)
		: Hex(id, type, nodes), omega(0.0) {}
	InfHex(const InfHex& other)
		: Hex(other), omega(other.omega) {}
	InfHex& operator=(const InfHex& other) {
		if (this != &other) {
			Hex::operator=(other);
			omega = other.omega;
		}
		return *this;
	}
	InfHex(InfHex&& other) noexcept
		: Hex(std::move(other)), omega(other.omega) {}
	InfHex& operator=(InfHex&& other) noexcept {
		if (this != &other) {
			Hex::operator=(std::move(other));
			omega = other.omega;
		}
		return *this;
	}
	virtual ~InfHex() = default;

	double omega = 0.0;
	
	double pole_x = 0.0, pole_y = 0.0, pole_z = 0.0;
	void set_pole(double px, double py, double pz) {
		pole_x = px; pole_y = py; pole_z = pz;
	}

	std::vector<std::complex<double>> FF(double ksi, double eta, double zeta) override;
	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;

protected:
	Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta) override;

private:
	double compute_A() const;
};
