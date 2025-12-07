// #define DISABLE_MATHMV_TESTS

#ifndef DISABLE_MATHMV_TESTS

#include "MathMVTests.h"
#include "Tri.h"
#include "Quad.h"

TEST_F(MathMVTest, BerlageFunction_ZeroTime) {
    double omega = 10.0;
    double A = 1.0;
    double time = 0.0;
    
    double result = berlage(omega, A, time);
    EXPECT_NEAR(result, 0.0, 1e-6);
}

TEST_F(MathMVTest, BerlageFunction_NegativeAmplitude) {
    double omega = 10.0;
    double A = -1.0;
    double time = 0.1;
    
    double result_neg = berlage(omega, A, time);
    double result_pos = berlage(omega, -A, time);
    
    EXPECT_NEAR(result_neg, -result_pos, tolerance);
}


TEST_F(GLLNodesTest, Order3_FourPoints) {
    std::vector<double> points, weights;
    compute_gll_nodes_weights(3, points, weights);
    
    ASSERT_EQ(points.size(), 4);
    ASSERT_EQ(weights.size(), 4);
    
    EXPECT_NEAR(points[0], -1.0, tolerance);
    EXPECT_NEAR(points[3], 1.0, tolerance);
    
    EXPECT_TRUE(verifyWeightsSum(weights));
    EXPECT_TRUE(verifyNodesSymmetry(points));
}

TEST_F(GLLNodesTest, Order4_FivePoints) {
    std::vector<double> points, weights;
    compute_gll_nodes_weights(4, points, weights);
    
    ASSERT_EQ(points.size(), 5);
    ASSERT_EQ(weights.size(), 5);
    
    EXPECT_NEAR(points[0], -1.0, tolerance);
    EXPECT_NEAR(points[4], 1.0, tolerance);
    
    EXPECT_NEAR(points[2], 0.0, tolerance);
    
    EXPECT_TRUE(verifyWeightsSum(weights));
    EXPECT_TRUE(verifyNodesSymmetry(points));
}

TEST_F(GLLNodesTest, AllWeightsPositive) {
    for (int order = 1; order <= 8; order++) {
        std::vector<double> points, weights;
        compute_gll_nodes_weights(order, points, weights);
        
        for (const auto& w : weights) {
            EXPECT_GT(w, 0.0) << "Order " << order << " has non-positive weight";
        }
    }
}

TEST_F(GLLNodesTest, PointsInRange) {
    for (int order = 1; order <= 8; order++) {
        std::vector<double> points, weights;
        compute_gll_nodes_weights(order, points, weights);
        
        for (const auto& p : points) {
            EXPECT_GE(p, -1.0) << "Order " << order << " has point less than -1";
            EXPECT_LE(p, 1.0) << "Order " << order << " has point greater than 1";
        }
    }
}

TEST_F(GLLNodesTest, PointsAreSorted) {
    for (int order = 1; order <= 8; order++) {
        std::vector<double> points, weights;
        compute_gll_nodes_weights(order, points, weights);
        
        for (size_t i = 1; i < points.size(); i++) {
            EXPECT_LT(points[i - 1], points[i]) 
                << "Order " << order << " points are not sorted at index " << i;
        }
    }
}

TEST_F(GaussRadauTest, Order4_FourPoints) {
    std::vector<double> points, weights;
    compute_gauss_radau_nodes_weights(4, points, weights);
    
    ASSERT_EQ(points.size(), 4);
    ASSERT_EQ(weights.size(), 4);
    
    EXPECT_NEAR(points[0], -1.0, tolerance);
}

TEST_F(GaussRadauTest, Order5_FivePoints) {
    std::vector<double> points, weights;
    compute_gauss_radau_nodes_weights(5, points, weights);
    
    ASSERT_EQ(points.size(), 5);
    ASSERT_EQ(weights.size(), 5);
    
    EXPECT_NEAR(points[0], -1.0, tolerance);
}

TEST_F(GaussRadauTest, FirstPointAlwaysMinusOne) {
    for (int n = 1; n <= 5; n++) {
        std::vector<double> points, weights;
        compute_gauss_radau_nodes_weights(n, points, weights);
        
        EXPECT_NEAR(points[0], -1.0, tolerance) 
            << "First point should be -1 for n=" << n;
    }
}

TEST_F(GaussRadauTest, AllWeightsPositive) {
    for (int n = 1; n <= 5; n++) {
        std::vector<double> points, weights;
        compute_gauss_radau_nodes_weights(n, points, weights);
        
        for (const auto& w : weights) {
            EXPECT_GT(w, 0.0) << "n=" << n << " has non-positive weight";
        }
    }
}

TEST_F(GaussRadauTest, PointsInRange) {
    for (int n = 1; n <= 5; n++) {
        std::vector<double> points, weights;
        compute_gauss_radau_nodes_weights(n, points, weights);
        
        for (const auto& p : points) {
            EXPECT_GE(p, -1.0) << "n=" << n << " has point less than -1";
            EXPECT_LE(p, 1.0) << "n=" << n << " has point greater than 1";
        }
    }
}

TEST_F(MathMVTest, BerlageFunction_ExactValues) {
    double omega = 10.0;
    double A = 1.0;
    
    double result_t0 = berlage(omega, A, 0.0);
    EXPECT_NEAR(result_t0, 0.0, 1e-10);
    
    double time = 0.01;
    double pi = 3.14159265358979323;
    double omega0 = 2 * pi * omega;
    double omega1 = omega0 / std::sqrt(3.0);
    double mult = A * std::pow(omega1, 2) * std::exp(-omega1 * time) / 4.0;
    double term1 = std::sin(omega0 * time) * (-1 * std::pow(time, 2) / omega1 + time / std::pow(omega1, 2) + 1 / std::pow(omega1, 3));
    double term2 = std::cos(omega0 * time) * (std::pow(time, 2) / omega1 + time / std::pow(omega1, 2));
    double expected = mult * (term1 - std::sqrt(3.0) * term2);
    
    double result = berlage(omega, A, time);
    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(MathMVTest, BerlageFunction_Scaling) {
    double omega = 10.0;
    double time = 0.1;
    
    double result_A1 = berlage(omega, 1.0, time);
    double result_A2 = berlage(omega, 2.0, time);
    
    EXPECT_NEAR(result_A2, 2.0 * result_A1, 1e-9);
}

TEST_F(MathMVTest, MinEdgeLength_ExactCalculation) {
    std::vector<int> nodes = { 1, 2, 3 };
    auto tri = std::make_shared<Tri>(1, ElemType::TRI, nodes);
    
    std::vector<double> x = { 0.0, 2.0, 0.0 };
    std::vector<double> y = { 0.0, 0.0, 1.0 };
    std::vector<double> z = { 0.0, 0.0, 0.0 };
    tri->set_coords(x, y, z);
    
    double min_len = min_edge_length(tri);
    
    EXPECT_NEAR(min_len, 1.0, 1e-9);
}

TEST_F(MathMVTest, MinEdgeLength_WithDifferentEdges) {
    std::vector<int> nodes = { 1, 2, 3, 4 };
    auto quad = std::make_shared<Quad>(1, ElemType::QUAD, nodes);
    
    std::vector<double> x = { 0.0, 1.0, 2.0, 0.0 };
    std::vector<double> y = { 0.0, 0.0, 2.0, 1.0 };
    std::vector<double> z = { 0.0, 0.0, 0.0, 0.0 };
    quad->set_coords(x, y, z);
    
    double min_len = min_edge_length(quad);
    
    EXPECT_NEAR(min_len, 1.0, 1e-9);
}

TEST_F(GLLNodesTest, GLLNodes_ExactValuesOrder2) {
    std::vector<double> points, weights;
    compute_gll_nodes_weights(2, points, weights);
    
    EXPECT_NEAR(points[0], -1.0, 1e-14);
    EXPECT_NEAR(points[1], 0.0, 1e-14);
    EXPECT_NEAR(points[2], 1.0, 1e-14);
    
    EXPECT_NEAR(weights[0], 1.0 / 3.0, 1e-14);
    EXPECT_NEAR(weights[1], 4.0 / 3.0, 1e-14);
    EXPECT_NEAR(weights[2], 1.0 / 3.0, 1e-14);
}

TEST_F(GLLNodesTest, GLLNodes_ExactValuesOrder3) {
    std::vector<double> points, weights;
    compute_gll_nodes_weights(3, points, weights);
    
    EXPECT_NEAR(points[0], -1.0, 1e-14);
    EXPECT_NEAR(points[1], -std::sqrt(1.0 / 5.0), 1e-14);
    EXPECT_NEAR(points[2], std::sqrt(1.0 / 5.0), 1e-14);
    EXPECT_NEAR(points[3], 1.0, 1e-14);
    
    double sum = 0.0;
    for (const auto& w : weights) {
        sum += w;
    }
    EXPECT_NEAR(sum, 2.0, 1e-14);
}

TEST_F(GaussRadauTest, GaussRadau_ExactValuesOrder2) {
    std::vector<double> points, weights;
    compute_gauss_radau_nodes_weights(2, points, weights);
    
    EXPECT_NEAR(points[0], -1.0, 1e-14);
    EXPECT_NEAR(points[1], 1.0 / 3.0, 1e-14);
    
    EXPECT_NEAR(weights[0], 0.5, 1e-14);
    EXPECT_NEAR(weights[1], 1.5, 1e-14);
}

TEST_F(GaussRadauTest, GaussRadau_ExactValuesOrder3) {
    std::vector<double> points, weights;
    compute_gauss_radau_nodes_weights(3, points, weights);
    
    double sqrt6 = std::sqrt(6.0);
    EXPECT_NEAR(points[0], -1.0, 1e-14);
    EXPECT_NEAR(points[1], (1.0 - sqrt6) / 5.0, 1e-14);
    EXPECT_NEAR(points[2], (1.0 + sqrt6) / 5.0, 1e-14);
    
    EXPECT_NEAR(weights[0], 2.0 / 9.0, 1e-14);
    EXPECT_NEAR(weights[1], (16.0 + sqrt6) / 18.0, 1e-14);
    EXPECT_NEAR(weights[2], (16.0 - sqrt6) / 18.0, 1e-14);
}

TEST_F(GaussRadauTest, GaussRadau_WeightsSum) {
    for (int n = 1; n <= 5; n++) {
        std::vector<double> points, weights;
        compute_gauss_radau_nodes_weights(n, points, weights);
        
        double sum = 0.0;
        for (const auto& w : weights) {
            sum += w;
        }
        EXPECT_NEAR(sum, 2.0, 1e-12) << "Sum of weights should be 2 for n=" << n;
    }
}

#endif
