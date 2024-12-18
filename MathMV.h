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
#include "log.h"

// MV - matrix vector

void printMatrix(Eigen::SparseMatrix<double> A);
void printVector(Eigen::SparseVector<double> B);
void productMV(Eigen::MatrixXd A, std::vector<double> b, std::vector<double>& c);
double line(std::vector<double> a, std::vector<double> b, std::vector<double> x);
void output_points(std::vector<double> start, std::vector<double> end, int count, std::vector<std::vector<double>>& points);

void KirshError(std::string filename, double P, double R);
void meshConvergence(std::string name1, std::string name2, std::string name3);

void field_name(std::string& filename, int type, int comp, int dim);
int output_fields(int type, int dim);
