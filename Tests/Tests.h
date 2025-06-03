#include <gtest/gtest.h>

#include "Node/Node.h"
#include "Element.h"
#include "Tri.h"
#include "Quad.h"

int Test();

class NodeTest : public ::testing::Test {
protected:
    void SetUp() override {
        coords1 = { 1.0, 2.0, 3.0 };
        coords2 = { 4.0, 5.0, 6.0 };
        node1 = std::make_unique<Node>(1, coords1);
        node2 = std::make_unique<Node>(2, coords2);
    }

    void TearDown() override {
        node1.reset();
        node2.reset();
    }

    std::array<double, 3> coords1;
    std::array<double, 3> coords2;
    std::unique_ptr<Node> node1;
    std::unique_ptr<Node> node2;
};

inline ::testing::AssertionResult ComplexNear(
    const std::complex<double>& a,
    const std::complex<double>& b,
    double tolerance = 1e-9) {
    if (std::abs(a.real() - b.real()) < tolerance &&
        std::abs(a.imag() - b.imag()) < tolerance) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
        << "Expected: (" << a.real() << ", " << a.imag() << ")\n"
        << "Actual: (" << b.real() << ", " << b.imag() << ")";
}

class TriTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<int> nodes = { 1, 2, 3 };
        tri = std::make_unique<Tri>(1, ElemType::TRI, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0 };
        tri->set_coords(x_coords, y_coords, z_coords);

        tri->set_constants(210e9, 0.3, 7850.0);

        expected_area = 0.5;
    }

    std::unique_ptr<Tri> tri;
    double expected_area;
};

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
        inf_quad = std::make_unique<infQuad>(1, ElemType::INFQUAD, nodes);

        std::vector<double> x_coords = { 0.0, 1.0, 1.0, 0.0 };
        std::vector<double> y_coords = { 0.0, 0.0, 1.0, 1.0 };
        std::vector<double> z_coords = { 0.0, 0.0, 0.0, 0.0 };
        inf_quad->set_coords(x_coords, y_coords, z_coords);

        inf_quad->set_constants(210e9, 0.3, 7850.0);
        inf_quad->is_dyn = true;
        inf_quad->omega = 100.0;
    }

    std::unique_ptr<infQuad> inf_quad;
};
