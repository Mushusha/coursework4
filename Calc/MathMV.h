#pragma once
#include <iostream>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"

#include "Enums.h"
#include "Element.h"
#include "Parser.h"
#include "Data.h"
#include "log1.h"

// MV - matrix vector

void print(Eigen::SparseMatrix<std::complex<double>>& A);
void print(Eigen::SparseVector<std::complex<double>>& B);
void print(Eigen::VectorX <std::complex<double>>& B);
double line(std::vector<double> a, std::vector<double> b, std::vector<double> x);

double berlage(double omega, double A, double time);

void compute_gll_nodes_weights(int order, std::vector<double>& points, std::vector<double>& weights);
