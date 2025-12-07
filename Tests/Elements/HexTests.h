#pragma once

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

#include "Hex.h"
#include "InfHex.h"
#include "Enums.h"

class HexTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4, 5, 6, 7, 8 };
        hex = std::make_unique<Hex>(1, ElemType::HEX, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
        hex->set_coords(x_coords, y_coords, z_coords);

        hex->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        hex->D = Eigen::MatrixXd::Zero(6, 6);
        hex->D(0, 0) = hex->D(1, 1) = hex->D(2, 2) = factor * (1.0 - nu);
        hex->D(0, 1) = hex->D(0, 2) = hex->D(1, 0) = hex->D(1, 2) = hex->D(2, 0) = hex->D(2, 1) = factor * nu;
        hex->D(3, 3) = hex->D(4, 4) = hex->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 1.0;
    }

    std::unique_ptr<Hex> hex;
    double expected_volume;
};

class HexNonUnitTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4, 5, 6, 7, 8 };
        hex = std::make_unique<Hex>(2, ElemType::HEX, nodes);

        std::vector<double> x_coords = { 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 3.0, 3.0, 0.0, 0.0, 3.0, 3.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0 };
        hex->set_coords(x_coords, y_coords, z_coords);

        hex->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        hex->D = Eigen::MatrixXd::Zero(6, 6);
        hex->D(0, 0) = hex->D(1, 1) = hex->D(2, 2) = factor * (1.0 - nu);
        hex->D(0, 1) = hex->D(0, 2) = hex->D(1, 0) = hex->D(1, 2) = hex->D(2, 0) = hex->D(2, 1) = factor * nu;
        hex->D(3, 3) = hex->D(4, 4) = hex->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;

        expected_volume = 24.0;
    }

    std::unique_ptr<Hex> hex;
    double expected_volume;
};

class InfHexTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4, 5, 6, 7, 8 };
        inf_hex = std::make_unique<InfHex>(1, ElemType::INFHEX, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
        inf_hex->set_coords(x_coords, y_coords, z_coords);

        inf_hex->set_constants(210e9, 0.3, 7850.0);
        inf_hex->is_dyn = false;
        inf_hex->omega = 100.0;
    }

    std::unique_ptr<InfHex> inf_hex;
};
