#pragma once

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "MathMV.h"

class MathMVTest : public ::testing::Test {
protected:
    void SetUp() override {
        tolerance = 1e-9;
    }

    double tolerance;
};

class GLLNodesTest : public ::testing::Test {
protected:
    void SetUp() override {
        tolerance = 1e-12;
    }

    double tolerance;

    bool verifyWeightsSum(const std::vector<double>& weights) {
        double sum = 0.0;
        for (const auto& w : weights) {
            sum += w;
        }
        return std::abs(sum - 2.0) < tolerance;
    }

    bool verifyNodesSymmetry(const std::vector<double>& points) {
        int n = points.size();
        for (int i = 0; i < n / 2; i++) {
            if (std::abs(points[i] + points[n - 1 - i]) > tolerance) {
                return false;
            }
        }
        return true;
    }
};

class GaussRadauTest : public ::testing::Test {
protected:
    void SetUp() override {
        tolerance = 1e-10;
    }

    double tolerance;
};
