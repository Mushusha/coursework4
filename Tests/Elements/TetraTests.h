#pragma once

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

#include "Tetra.h"
#include "Enums.h"

class TetraTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4 };
        tetra = std::make_unique<Tetra>(1, ElemType::TETRA, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 0.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 0.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 1.0 };
        tetra->set_coords(x_coords, y_coords, z_coords);

        tetra->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        tetra->D = Eigen::MatrixXd::Zero(6, 6);
        tetra->D(0, 0) = tetra->D(1, 1) = tetra->D(2, 2) = factor * (1.0 - nu);
        tetra->D(0, 1) = tetra->D(0, 2) = tetra->D(1, 0) = tetra->D(1, 2) = tetra->D(2, 0) = tetra->D(2, 1) = factor * nu;
        tetra->D(3, 3) = tetra->D(4, 4) = tetra->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 1.0 / 6.0;
    }

    std::unique_ptr<Tetra> tetra;
    double expected_volume;
};

class TetraScaledTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4 };
        tetra = std::make_unique<Tetra>(2, ElemType::TETRA, nodes);

        std::vector<double> x_coords = { 0.0, 2.0, 0.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 2.0, 0.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 2.0 };
        tetra->set_coords(x_coords, y_coords, z_coords);

        tetra->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        tetra->D = Eigen::MatrixXd::Zero(6, 6);
        tetra->D(0, 0) = tetra->D(1, 1) = tetra->D(2, 2) = factor * (1.0 - nu);
        tetra->D(0, 1) = tetra->D(0, 2) = tetra->D(1, 0) = tetra->D(1, 2) = tetra->D(2, 0) = tetra->D(2, 1) = factor * nu;
        tetra->D(3, 3) = tetra->D(4, 4) = tetra->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 4.0 / 3.0;
    }

    std::unique_ptr<Tetra> tetra;
    double expected_volume;
};
