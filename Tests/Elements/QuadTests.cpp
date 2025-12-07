// #define DISABLE_QUAD_TESTS

#ifndef DISABLE_QUAD_TESTS

#include "QuadTests.h"
#include <Eigen/Eigenvalues>

TEST_F(QuadTest, FF_ExactValues) {

    auto ff_node1 = quad->FF(-1.0, -1.0, 0.0);
    EXPECT_NEAR(std::abs(ff_node1[0]), 1.0, 1e-9) << "FF[0] should be equal to 1.0";
    EXPECT_NEAR(std::abs(ff_node1[1]), 0.0, 1e-9) << "FF[1] should be equal to 0.0";
    EXPECT_NEAR(std::abs(ff_node1[2]), 0.0, 1e-9) << "FF[2] should be equal to 0.0";
    EXPECT_NEAR(std::abs(ff_node1[3]), 0.0, 1e-9) << "FF[3] should be equal to 0.0";

    auto ff_center = quad->FF(0.0, 0.0, 0.0);
    const double expected_center = 0.25;
    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(std::abs(ff_center[i]), expected_center, 1e-9)
            << "FF[" << i << "] at center should be equal to " << expected_center;
    }
}

TEST_F(QuadTest, B_ExactAllValuesAtCenter) {

    Eigen::MatrixXcd B = quad->B(0.0, 0.0);
    
    const double expected_dN_X0 = -0.5;
    const double expected_dN_Y0 = -0.5;
    
    EXPECT_NEAR(B(0, 0).real(), expected_dN_X0, 1e-8) << "B(0,0) should be equal to " << expected_dN_X0;
    EXPECT_NEAR(B(1, 1).real(), expected_dN_Y0, 1e-8) << "B(1,1) should be equal to " << expected_dN_Y0;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 8; j++) {
            EXPECT_TRUE(std::isfinite(B(i, j).real()))
                << "B(" << i << "," << j << ") should be finite";
            EXPECT_NEAR(B(i, j).imag(), 0.0, 1e-12)
                << "B(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    EXPECT_EQ(B.rows(), 3) << "B should have 3 rows";
    EXPECT_EQ(B.cols(), 8) << "B should have 8 columns";
}

TEST_F(QuadTest, localK_ExactAllMatrixValues) {

    double E = 210e9;
    double nu = 0.3;
    quad->set_constants(E, nu, 7850.0);
    
    double factor = E / (1 - nu * nu);
    quad->D = Eigen::MatrixXd::Zero(3, 3);
    quad->D(0, 0) = quad->D(1, 1) = factor;
    quad->D(0, 1) = quad->D(1, 0) = factor * nu;
    quad->D(2, 2) = factor * (1 - nu) / 2.0;
    
    Eigen::MatrixXcd K = quad->localK();
    
    const double K_values[8][8] = {
        {230769230769.23077, 69230769230.76923, -115384615384.61539, -34615384615.38462, -115384615384.61539, -34615384615.38462, 0.0, 0.0},
        {69230769230.76923, 230769230769.23077, -34615384615.38462, -115384615384.61539, -34615384615.38462, -115384615384.61539, 0.0, 0.0},
        {-115384615384.61539, -34615384615.38462, 230769230769.23077, 69230769230.76923, 0.0, 0.0, -115384615384.61539, -34615384615.38462},
        {-34615384615.38462, -115384615384.61539, 69230769230.76923, 230769230769.23077, 0.0, 0.0, -34615384615.38462, -115384615384.61539},
        {-115384615384.61539, -34615384615.38462, 0.0, 0.0, 230769230769.23077, 69230769230.76923, -115384615384.61539, -34615384615.38462},
        {-34615384615.38462, -115384615384.61539, 0.0, 0.0, 69230769230.76923, 230769230769.23077, -34615384615.38462, -115384615384.61539},
        {0.0, 0.0, -115384615384.61539, -34615384615.38462, -115384615384.61539, -34615384615.38462, 230769230769.23077, 69230769230.76923},
        {0.0, 0.0, -34615384615.38462, -115384615384.61539, -34615384615.38462, -115384615384.61539, 69230769230.76923, 230769230769.23077}
    };
    
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            const double expected = K_values[i][j];
            EXPECT_NEAR(K(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "K(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(K(i, j).imag(), 0.0, 1e-12)
                << "K(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            EXPECT_NEAR(K(i, j).real(), K(j, i).real(), std::abs(K(i, j).real()) * 1e-10)
                << "K is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 8; i++) {
        EXPECT_GT(K(i, i).real(), 0.0)
            << "Diagonal element K(" << i << "," << i << ") should be positive";
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(K.real());
    for (int i = 0; i < 8; i++) {
        EXPECT_GT(solver.eigenvalues()(i), 0.0)
            << "Eigenvalue " << i << " should be positive";
    }
}

TEST_F(QuadTest, localF_ExactAllValues) {

    double pressure_value = 1000.0;
    quad->set_load(PRESSURE, 0, {pressure_value, 0, 0, 0, 0, 0});
    
    std::vector<double> F = quad->localF();
    
    const double F_expected[8] = {
        0.0, 500.0, 0.0, 500.0, 0.0, 0.0, 0.0, 0.0
    };
    
    EXPECT_EQ(F.size(), 8) << "Vector size F should be 8";
    for (int i = 0; i < 8; i++) {
        EXPECT_NEAR(F[i], F_expected[i], std::abs(F_expected[i]) * 1e-6)
            << "F[" << i << "] should be equal to " << F_expected[i];
        EXPECT_TRUE(std::isfinite(F[i])) << "F[" << i << "] should be finite";
    }
}

TEST_F(QuadTest, localC_ExactAllMatrixValues) {

    Eigen::MatrixXd C = quad->localC();
    
    const double C_values[4][4] = {
        {0.111111111111111, 0.0555555555555556, 0.0277777777777778, 0.0555555555555555},
        {0.0555555555555556, 0.111111111111111, 0.0555555555555556, 0.0277777777777778},
        {0.0277777777777778, 0.0555555555555556, 0.111111111111111, 0.0555555555555555},
        {0.0555555555555555, 0.0277777777777778, 0.0555555555555555, 0.111111111111111}
    };
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            const double expected = C_values[i][j];
            EXPECT_NEAR(C(i, j), expected, std::abs(expected) * 1e-10)
                << "C(" << i << "," << j << ") should be equal to " << expected;
        }
    }
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_NEAR(C(i, j), C(j, i), 1e-12)
                << "C is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_GT(C(i, j), 0.0)
                << "C(" << i << "," << j << ") should be positive";
        }
    }
}

TEST_F(QuadTest, localR_ExactAllValues) {

    std::vector<double> values = { 1.0, 1.0, 1.0, 1.0 };
    std::vector<double> R = quad->localR(values);
    
    const double expected = 0.25;
    
    EXPECT_EQ(R.size(), 4) << "Vector size R should be 4";
    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(R[i], expected, 1e-9)
            << "R[" << i << "] should be equal to " << expected;
    }
    
    for (int i = 0; i < 4; i++) {
        EXPECT_GT(R[i], 0.0) << "R[" << i << "] should be positive";
    }
}

TEST_F(QuadTest, localM_ExactAllMatrixValues) {

    Eigen::MatrixXcd M = quad->localM();
    
    const double expected_diag = 1962.5;
    
    for (int i = 0; i < 8; i += 2) {
        EXPECT_NEAR(M(i, i).real(), expected_diag, 1e-6)
            << "M(" << i << "," << i << ") should be equal to " << expected_diag;
    }
    for (int i = 1; i < 8; i += 2) {
        EXPECT_NEAR(M(i, i).real(), expected_diag, 1e-6)
            << "M(" << i << "," << i << ") should be equal to " << expected_diag;
    }
    
    const double M_values[8][8] = {
        {872.222222222222, 0.0, 436.111111111111, 0.0, 218.055555555555, 0.0, 436.111111111111, 0.0},
        {0.0, 872.222222222222, 0.0, 436.111111111111, 0.0, 218.055555555555, 0.0, 436.111111111111},
        {436.111111111111, 0.0, 872.222222222222, 0.0, 436.111111111111, 0.0, 218.055555555555, 0.0},
        {0.0, 436.111111111111, 0.0, 872.222222222222, 0.0, 436.111111111111, 0.0, 218.055555555555},
        {218.055555555555, 0.0, 436.111111111111, 0.0, 872.222222222222, 0.0, 436.111111111111, 0.0},
        {0.0, 218.055555555555, 0.0, 436.111111111111, 0.0, 872.222222222222, 0.0, 436.111111111111},
        {436.111111111111, 0.0, 218.055555555555, 0.0, 436.111111111111, 0.0, 872.222222222222, 0.0},
        {0.0, 436.111111111111, 0.0, 218.055555555555, 0.0, 436.111111111111, 0.0, 872.222222222222}
    };
    
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            const double expected = M_values[i][j];
            EXPECT_NEAR(M(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "M(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(M(i, j).imag(), 0.0, 1e-12)
                << "M(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            EXPECT_NEAR(M(i, j).real(), M(j, i).real(), std::abs(M(i, j).real()) * 1e-10)
                << "M is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 8; i++) {
        EXPECT_GT(M(i, i).real(), 0.0)
            << "Diagonal element M(" << i << "," << i << ") should be positive";
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M.real());
    for (int i = 0; i < 8; i++) {
        EXPECT_GT(solver.eigenvalues()(i), 0.0)
            << "Eigenvalue " << i << " should be positive";
    }
}

TEST_F(QuadTest, Volume_ExactAreaCalculation) {

    double area = quad->Volume();
    const double expected = 1.0;
    EXPECT_NEAR(area, expected, 1e-9) << "Area should be equal to " << expected;
}

#endif
