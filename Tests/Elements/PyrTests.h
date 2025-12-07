#pragma once

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

#include "Pyr.h"
#include "Enums.h"

class PyrTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4, 5 };
        pyr = std::make_unique<Pyr>(1, ElemType::PYR, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 1.0, 0.0, 0.5 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 1.0, 0.5 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0, 1.0 };
        pyr->set_coords(x_coords, y_coords, z_coords);

        pyr->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        pyr->D = Eigen::MatrixXd::Zero(6, 6);
        pyr->D(0, 0) = pyr->D(1, 1) = pyr->D(2, 2) = factor * (1.0 - nu);
        pyr->D(0, 1) = pyr->D(0, 2) = pyr->D(1, 0) = pyr->D(1, 2) = pyr->D(2, 0) = pyr->D(2, 1) = factor * nu;
        pyr->D(3, 3) = pyr->D(4, 4) = pyr->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 1.0 / 3.0;
    }

    std::unique_ptr<Pyr> pyr;
    double expected_volume;
};

class PyrScaledTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4, 5 };
        pyr = std::make_unique<Pyr>(2, ElemType::PYR, nodes);

        std::vector<double> x_coords = { 0.0, 2.0, 2.0, 0.0, 1.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 3.0, 3.0, 1.5 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0, 4.0 };
        pyr->set_coords(x_coords, y_coords, z_coords);

        pyr->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        pyr->D = Eigen::MatrixXd::Zero(6, 6);
        pyr->D(0, 0) = pyr->D(1, 1) = pyr->D(2, 2) = factor * (1.0 - nu);
        pyr->D(0, 1) = pyr->D(0, 2) = pyr->D(1, 0) = pyr->D(1, 2) = pyr->D(2, 0) = pyr->D(2, 1) = factor * nu;
        pyr->D(3, 3) = pyr->D(4, 4) = pyr->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 8.0;
    }

    std::unique_ptr<Pyr> pyr;
    double expected_volume;
};
