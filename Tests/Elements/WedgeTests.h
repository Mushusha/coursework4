#pragma once

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

#include "Wedge.h"
#include "Enums.h"

class WedgeTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4, 5, 6 };
        wedge = std::make_unique<Wedge>(1, ElemType::WEDGE, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
        wedge->set_coords(x_coords, y_coords, z_coords);

        wedge->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        wedge->D = Eigen::MatrixXd::Zero(6, 6);
        wedge->D(0, 0) = wedge->D(1, 1) = wedge->D(2, 2) = factor * (1.0 - nu);
        wedge->D(0, 1) = wedge->D(0, 2) = wedge->D(1, 0) = wedge->D(1, 2) = wedge->D(2, 0) = wedge->D(2, 1) = factor * nu;
        wedge->D(3, 3) = wedge->D(4, 4) = wedge->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 0.5;
    }

    std::unique_ptr<Wedge> wedge;
    double expected_volume;
};

class WedgeScaledTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4, 5, 6 };
        wedge = std::make_unique<Wedge>(2, ElemType::WEDGE, nodes);

        std::vector<double> x_coords = { 0.0, 2.0, 0.0, 0.0, 2.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 2.0, 0.0, 0.0, 2.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 3.0, 3.0, 3.0 };
        wedge->set_coords(x_coords, y_coords, z_coords);

        wedge->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        wedge->D = Eigen::MatrixXd::Zero(6, 6);
        wedge->D(0, 0) = wedge->D(1, 1) = wedge->D(2, 2) = factor * (1.0 - nu);
        wedge->D(0, 1) = wedge->D(0, 2) = wedge->D(1, 0) = wedge->D(1, 2) = wedge->D(2, 0) = wedge->D(2, 1) = factor * nu;
        wedge->D(3, 3) = wedge->D(4, 4) = wedge->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 6.0;
    }

    std::unique_ptr<Wedge> wedge;
    double expected_volume;
};
