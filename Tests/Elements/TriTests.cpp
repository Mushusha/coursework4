// #define DISABLE_TRI_TESTS

#ifndef DISABLE_TRI_TESTS

#include "TriTests.h"
#include <Eigen/Eigenvalues>

TEST_F(TriTest, FF_ExactValues) {

    auto ff_node1 = tri->FF(0.3, 0.0, 0.0); 
    EXPECT_NEAR(std::abs(ff_node1[0]), 0.7, 1e-9) << "FF[0] should be equal to 0.7";
    EXPECT_NEAR(std::abs(ff_node1[1]), 0.3, 1e-9) << "FF[1] should be equal to 0.3";
    EXPECT_NEAR(std::abs(ff_node1[2]), 0.0, 1e-9) << "FF[2] should be equal to 0.0";

    auto ff_node2 = tri->FF(0.0, 0.3, 0.0); 
    EXPECT_NEAR(std::abs(ff_node2[0]), 0.7, 1e-9) << "FF[0] should be equal to 0.7";
    EXPECT_NEAR(std::abs(ff_node2[1]), 0.0, 1e-9) << "FF[1] should be equal to 0.0";
    EXPECT_NEAR(std::abs(ff_node2[2]), 0.3, 1e-9) << "FF[2] should be equal to 0.3";

    auto ff_center = tri->FF(1.0 / 3.0, 1.0 / 3.0, 0.0);
    const double expected_center = 0.3333333333333333;
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(std::abs(ff_center[i]), expected_center, 1e-9)
            << "FF[" << i << "] at centroid should be equal to " << expected_center;
    }
}

TEST_F(TriTest, B_ExactAllValues) {

    Eigen::MatrixXcd B = tri->B();
    
    const double expected_B00 = -1.0;
    const double expected_B11 = -1.0;
    const double expected_B20 = -1.0;
    const double expected_B21 = -1.0;
    
    EXPECT_NEAR(B(0, 0).real(), expected_B00, 1e-9) << "B(0,0) should be equal to " << expected_B00;
    EXPECT_NEAR(B(1, 1).real(), expected_B11, 1e-9) << "B(1,1) should be equal to " << expected_B11;
    EXPECT_NEAR(B(2, 0).real(), expected_B20, 1e-9) << "B(2,0) should be equal to " << expected_B20;
    EXPECT_NEAR(B(2, 1).real(), expected_B21, 1e-9) << "B(2,1) should be equal to " << expected_B21;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            EXPECT_TRUE(std::isfinite(B(i, j).real()))
                << "B(" << i << "," << j << ") should be finite";
            EXPECT_NEAR(B(i, j).imag(), 0.0, 1e-12)
                << "B(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    EXPECT_EQ(B.rows(), 3) << "B should have 3 rows";
    EXPECT_EQ(B.cols(), 6) << "B should have 6 columns";
}

TEST_F(TriTest, localK_ExactAllMatrixValues) {

    double E = 210e9;
    double nu = 0.3;
    tri->set_constants(E, nu, 7850.0);
    
    double factor = E / (1 - nu * nu);
    tri->D = Eigen::MatrixXd::Zero(3, 3);
    tri->D(0, 0) = tri->D(1, 1) = factor;
    tri->D(0, 1) = tri->D(1, 0) = factor * nu;
    tri->D(2, 2) = factor * (1 - nu) / 2.0;
    
    Eigen::MatrixXcd K = tri->localK();
    
    const double K_values[6][6] = {
        {230769230769.23077, 69230769230.76923, -230769230769.23077, 0.0, 0.0, -69230769230.76923},
        {69230769230.76923, 230769230769.23077, -69230769230.76923, 0.0, -230769230769.23077, 0.0},
        {-230769230769.23077, -69230769230.76923, 230769230769.23077, 0.0, 0.0, 69230769230.76923},
        {0.0, 0.0, 0.0, 80769230769.23077, 69230769230.76923, -80769230769.23077},
        {0.0, -230769230769.23077, 0.0, 69230769230.76923, 230769230769.23077, -69230769230.76923},
        {-69230769230.76923, 0.0, 69230769230.76923, -80769230769.23077, -69230769230.76923, 80769230769.23077}
    };
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            const double expected = K_values[i][j];
            EXPECT_NEAR(K(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "K(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(K(i, j).imag(), 0.0, 1e-12)
                << "K(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            EXPECT_NEAR(K(i, j).real(), K(j, i).real(), std::abs(K(i, j).real()) * 1e-10)
                << "K is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 6; i++) {
        EXPECT_GT(K(i, i).real(), 0.0)
            << "Diagonal element K(" << i << "," << i << ") should be positive";
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(K.real());
    for (int i = 0; i < 6; i++) {
        EXPECT_GT(solver.eigenvalues()(i), 0.0)
            << "Eigenvalue " << i << " should be positive";
    }
}

TEST_F(TriTest, localF_ExactAllValues) {

    double pressure_value = 1000.0;
    tri->set_load(PRESSURE, 0, {pressure_value, 0, 0, 0, 0, 0});
    
    std::vector<double> F = tri->localF();
    
    const double F_expected[6] = {
        0.0, 500.0, 0.0, 500.0, 0.0, 0.0
    };
    
    EXPECT_EQ(F.size(), 6) << "Vector size F should be 6";
    for (int i = 0; i < 6; i++) {
        EXPECT_NEAR(F[i], F_expected[i], std::abs(F_expected[i]) * 1e-6)
            << "F[" << i << "] should be equal to " << F_expected[i];
        EXPECT_TRUE(std::isfinite(F[i])) << "F[" << i << "] should be finite";
    }
}

TEST_F(TriTest, localC_ExactAllMatrixValues) {

    Eigen::MatrixXd C = tri->localC();
    
    const double C_values[3][3] = {
        {0.0833333333333333, 0.0416666666666667, 0.0416666666666667},
        {0.0416666666666667, 0.0833333333333333, 0.0416666666666667},
        {0.0416666666666667, 0.0416666666666667, 0.0833333333333333}
    };
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            const double expected = C_values[i][j];
            EXPECT_NEAR(C(i, j), expected, std::abs(expected) * 1e-9)
                << "C(" << i << "," << j << ") should be equal to " << expected;
        }
    }
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_NEAR(C(i, j), C(j, i), 1e-12)
                << "C is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_GT(C(i, j), 0.0)
                << "C(" << i << "," << j << ") should be positive";
        }
    }
}

TEST_F(TriTest, localR_ExactAllValues) {

    std::vector<double> values = { 1.0, 1.0, 1.0 };
    std::vector<double> R = tri->localR(values);
    
    const double expected = 0.16666666666666666;
    
    EXPECT_EQ(R.size(), 3) << "Vector size R should be 3";
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(R[i], expected, 1e-9)
            << "R[" << i << "] should be equal to " << expected;
    }
    
    for (int i = 0; i < 3; i++) {
        EXPECT_GT(R[i], 0.0) << "R[" << i << "] should be positive";
    }
}

TEST_F(TriTest, localM_ExactAllMatrixValues) {

    tri->set_constants(210e9, 0.3, 7850.0);
    
    Eigen::MatrixXcd M = tri->localM();
    
    const double M_values[6][6] = {
        {654.166666666667, 0.0, 327.083333333333, 0.0, 327.083333333333, 0.0},
        {0.0, 654.166666666667, 0.0, 327.083333333333, 0.0, 327.083333333333},
        {327.083333333333, 0.0, 654.166666666667, 0.0, 327.083333333333, 0.0},
        {0.0, 327.083333333333, 0.0, 654.166666666667, 0.0, 327.083333333333},
        {327.083333333333, 0.0, 327.083333333333, 0.0, 654.166666666667, 0.0},
        {0.0, 327.083333333333, 0.0, 327.083333333333, 0.0, 654.166666666667}
    };
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            const double expected = M_values[i][j];
            EXPECT_NEAR(M(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "M(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(M(i, j).imag(), 0.0, 1e-12)
                << "M(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            EXPECT_NEAR(M(i, j).real(), M(j, i).real(), std::abs(M(i, j).real()) * 1e-10)
                << "M is not symmetric at position (" << i << "," << j << ")";
            EXPECT_NEAR(M(i, j).imag(), 0.0, 1e-12)
                << "M has non-zero imaginary part at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 6; i++) {
        EXPECT_GT(M(i, i).real(), 0.0)
            << "Diagonal element M(" << i << "," << i << ") should be positive";
    }
}

TEST_F(TriTest, Volume_ExactAreaCalculation) {

    double area = tri->Volume();
    const double expected = 0.5;
    EXPECT_NEAR(area, expected, 1e-9) << "Area should be equal to " << expected;
    
    std::vector<double> x_coords = { 0.0, 2.0, 0.0 };
    std::vector<double> y_coords = { 0.0, 0.0, 2.0 };
    tri->set_coords(x_coords, y_coords, { 0.0, 0.0, 0.0 });
    
    area = tri->Volume();
    const double expected_scaled = 2.0;
    EXPECT_NEAR(area, expected_scaled, 1e-9) 
        << "Scaled triangle area should be equal to " << expected_scaled;
}

#endif
