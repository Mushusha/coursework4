#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/Core"
#include "unsupported/Eigen/IterativeSolvers"

#include "Calc/Data.h"
#include "log1.h"


class Smoothing {
public:
	Smoothing(Data& data, Tensor::ResType type);
    Smoothing(const Smoothing& other);
    Smoothing& operator=(const Smoothing& other);
    Smoothing(Smoothing&& other) noexcept;
	Smoothing& operator=(Smoothing&& other) noexcept;
    virtual ~Smoothing() = default;

	void solve();
	//std::vector<double> line_start;
	//std::vector<double> line_end;
	//int points_count;
	//Smoothing(Data& data, std::vector<double> start,
	//	std::vector<double> end, int count);


private:
	Smoothing() = delete;

	Data calc_data;
	Tensor::ResType type;

	Eigen::SparseMatrix <double> C;
	Eigen::SparseVector <double> R;

	void fillGlobalC();
	void fillGlobalR(int type, int comp);
};