#pragma once

#include <gtest/gtest.h>
#include <memory>
#include <complex>

#include "Tri.h"
#include "Enums.h"

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
        << "Got: (" << b.real() << ", " << b.imag() << ")";
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
