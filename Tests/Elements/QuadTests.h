#pragma once

#include <gtest/gtest.h>
#include <memory>
#include <complex>

#include "Quad.h"
#include "InfQuad.h"
#include "Enums.h"

class QuadTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4 };
        quad = std::make_unique<Quad>(1, ElemType::QUAD, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 1.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 1.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0 };
        quad->set_coords(x_coords, y_coords, z_coords);

        quad->set_constants(210e9, 0.3, 7850.0);

        double E = 210e9;
        double nu = 0.3;
        quad->D = Eigen::MatrixXd::Zero(3, 3);
        quad->D << 1.0, nu, 0.0,
            nu, 1.0, 0.0,
            0.0, 0.0, (1.0 - nu) / 2.0;
        quad->D *= E / (1.0 - nu * nu);
    }

    std::unique_ptr<Quad> quad;
};

class infQuadTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3, 4 };
        inf_quad = std::make_unique<InfQuad>(1, ElemType::INFQUAD, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 1.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 1.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0 };
        inf_quad->set_coords(x_coords, y_coords, z_coords);

        inf_quad->set_constants(210e9, 0.3, 7850.0);
        inf_quad->omega = 100.0;
    }

    std::unique_ptr<InfQuad> inf_quad;
};
