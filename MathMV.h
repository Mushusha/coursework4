#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"

// MV - matrix vector

void writeMatrix(Eigen::SparseMatrix<double> A);
void writeVector(Eigen::SparseVector<double> B);
void productMV(Eigen::MatrixXd A, std::vector<double> b, std::vector<double>& c);
