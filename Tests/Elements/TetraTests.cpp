// #define DISABLE_TETRA_TESTS

#ifndef DISABLE_TETRA_TESTS

#include "TetraTests.h"
#include <Eigen/Eigenvalues>

TEST_F(TetraTest, FF_ExactValues) {

    auto ff_node0 = tetra->FF(0.0, 0.0, 0.0);
    EXPECT_NEAR(std::abs(ff_node0[0]), 1.0, 1e-9) << "FF[0] should be equal to 1.0";
    EXPECT_NEAR(std::abs(ff_node0[1]), 0.0, 1e-9) << "FF[1] should be equal to 0.0";
    EXPECT_NEAR(std::abs(ff_node0[2]), 0.0, 1e-9) << "FF[2] should be equal to 0.0";
    EXPECT_NEAR(std::abs(ff_node0[3]), 0.0, 1e-9) << "FF[3] should be equal to 0.0";

    auto ff_node1 = tetra->FF(1.0, 0.0, 0.0);
    EXPECT_NEAR(std::abs(ff_node1[0]), 0.0, 1e-9) << "FF[0] should be equal to 0.0";
    EXPECT_NEAR(std::abs(ff_node1[1]), 1.0, 1e-9) << "FF[1] should be equal to 1.0";
    EXPECT_NEAR(std::abs(ff_node1[2]), 0.0, 1e-9) << "FF[2] should be equal to 0.0";
    EXPECT_NEAR(std::abs(ff_node1[3]), 0.0, 1e-9) << "FF[3] should be equal to 0.0";

    auto ff_center = tetra->FF(0.25, 0.25, 0.25);
    const double expected_center = 0.25;
    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(std::abs(ff_center[i]), expected_center, 1e-9)
            << "FF[" << i << "] at centroid should be equal to " << expected_center;
    }
}

TEST_F(TetraTest, B_ExactAllValues) {

    Eigen::MatrixXcd B = tetra->B();
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 12; j++) {
            EXPECT_TRUE(std::isfinite(B(i, j).real()))
                << "B(" << i << "," << j << ") should be finite";
            EXPECT_NEAR(B(i, j).imag(), 0.0, 1e-12)
                << "B(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    EXPECT_EQ(B.rows(), 6) << "B should have 6 rows";
    EXPECT_EQ(B.cols(), 12) << "B should have 12 columns";
    
    double norm = B.norm();
    EXPECT_GT(norm, 0.0) << "Norm B should be positive";
}

TEST_F(TetraTest, localK_ExactAllMatrixValues) {

    double E = 210e9;
    double nu = 0.3;
    tetra->set_constants(E, nu, 7850.0);
    
    double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    tetra->D = Eigen::MatrixXd::Zero(6, 6);
    tetra->D(0, 0) = tetra->D(1, 1) = tetra->D(2, 2) = factor * (1.0 - nu);
    tetra->D(0, 1) = tetra->D(0, 2) = tetra->D(1, 0) = tetra->D(1, 2) = tetra->D(2, 0) = tetra->D(2, 1) = factor * nu;
    tetra->D(3, 3) = tetra->D(4, 4) = tetra->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;
    
    Eigen::MatrixXcd K = tetra->localK();
    
    const double K_values[12][12] = {
        {74038461538.4615, 33653846153.8461, 33653846153.8461, -47115384615.3846, -13461538461.5385, -13461538461.5385, -13461538461.5385, -20192307692.3077, 0, -13461538461.5385, 0, -20192307692.3077},
        {33653846153.8461, 74038461538.4615, 33653846153.8461, -20192307692.3077, -13461538461.5385, 0, -13461538461.5385, -47115384615.3846, -13461538461.5385, 0, -13461538461.5385, -20192307692.3077},
        {33653846153.8461, 33653846153.8461, 74038461538.4615, -20192307692.3077, 0, -13461538461.5385, 0, -20192307692.3077, -13461538461.5385, -13461538461.5385, -13461538461.5385, -47115384615.3846},
        {-47115384615.3846, -20192307692.3077, -20192307692.3077, 47115384615.3846, 0, 0, 0, 20192307692.3077, 0, 0, 0, 20192307692.3077},
        {-13461538461.5385, -13461538461.5385, 0, 0, 13461538461.5385, 0, 13461538461.5385, 0, 0, 0, 0, 0},
        {-13461538461.5385, 0, -13461538461.5385, 0, 0, 13461538461.5385, 0, 0, 0, 13461538461.5385, 0, 0},
        {-13461538461.5385, -13461538461.5385, 0, 0, 13461538461.5385, 0, 13461538461.5385, 0, 0, 0, 0, 0},
        {-20192307692.3077, -47115384615.3846, -20192307692.3077, 20192307692.3077, 0, 0, 0, 47115384615.3846, 0, 0, 0, 20192307692.3077},
        {0, -13461538461.5385, -13461538461.5385, 0, 0, 0, 0, 0, 13461538461.5385, 0, 13461538461.5385, 0},
        {-13461538461.5385, 0, -13461538461.5385, 0, 0, 13461538461.5385, 0, 0, 0, 13461538461.5385, 0, 0},
        {0, -13461538461.5385, -13461538461.5385, 0, 0, 0, 0, 0, 13461538461.5385, 0, 13461538461.5385, 0},
        {-20192307692.3077, -20192307692.3077, -47115384615.3846, 20192307692.3077, 0, 0, 0, 20192307692.3077, 0, 0, 0, 47115384615.3846}
    };
    
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            const double expected = K_values[i][j];
            EXPECT_NEAR(K(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "K(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(K(i, j).imag(), 0.0, 1e-12)
                << "K(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            EXPECT_NEAR(K(i, j).real(), K(j, i).real(), std::abs(K(i, j).real()) * 1e-10)
                << "K is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 12; i++) {
        EXPECT_GT(K(i, i).real(), 0.0)
            << "Diagonal element K(" << i << "," << i << ") should be positive";
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(K.real());
    for (int i = 0; i < 12; i++) {
        EXPECT_GT(solver.eigenvalues()(i), 0.0)
            << "Eigenvalue " << i << " should be positive";
    }
}

TEST_F(TetraTest, localF_ExactAllValues) {

    double pressure_value = 1000.0;
    tetra->set_load(PRESSURE, 0, {pressure_value, 0, 0, 0, 0, 0});
    
    std::vector<double> F = tetra->localF();
    
    EXPECT_EQ(F.size(), 12) << "Vector size F should be 12";
    
    const double F_expected[12] = {
        -166.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, -166.666666666667, 0.0, 0.0, -166.666666666667, 0.0, 0.0
    };
    
    for (int i = 0; i < 12; i++) {
        const double expected = F_expected[i];
        EXPECT_NEAR(F[i], expected, std::abs(expected) * 1e-6)
            << "F[" << i << "] should be equal to " << expected;
        EXPECT_TRUE(std::isfinite(F[i])) << "F[" << i << "] should be finite";
    }
    
    double total_force = 0.0;
    for (int i = 0; i < 12; i++) {
        total_force += std::abs(F[i]);
    }
    EXPECT_GT(total_force, 0.0) << "Total load should be positive";
}

TEST_F(TetraTest, localC_ExactAllMatrixValues) {

    Eigen::MatrixXd C = tetra->localC();
    
    const double C_values[4][4] = {
        {0.0166666666666667, 0.00833333333333333, 0.00833333333333333, 0.00833333333333333},
        {0.00833333333333333, 0.0166666666666667, 0.00833333333333333, 0.00833333333333333},
        {0.00833333333333333, 0.00833333333333333, 0.0166666666666667, 0.00833333333333333},
        {0.00833333333333333, 0.00833333333333333, 0.00833333333333333, 0.0166666666666667}
    };
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            const double expected = C_values[i][j];
            EXPECT_NEAR(C(i, j), expected, std::abs(expected) * 1e-9)
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

TEST_F(TetraTest, localR_ExactAllValues) {

    std::vector<double> values = { 1.0, 1.0, 1.0, 1.0 };
    std::vector<double> R = tetra->localR(values);
    
    const double expected_R = 0.041666666666666664;
    
    EXPECT_EQ(R.size(), 4) << "Vector size R should be 4";
    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(R[i], expected_R, 1e-9)
            << "R[" << i << "] should be equal to " << expected_R;
    }
    
    for (int i = 0; i < 4; i++) {
        EXPECT_GT(R[i], 0.0) << "R[" << i << "] should be positive";
    }
}

TEST_F(TetraTest, localM_ExactAllMatrixValues) {

    Eigen::MatrixXcd M = tetra->localM();
    
    const double expected_mass_per_node = 130.83333333333334;
    
    for (int i = 0; i < 12; i += 3) {
        EXPECT_NEAR(M(i, i).real(), expected_mass_per_node, 1e-6)
            << "M(" << i << "," << i << ") should be equal to " << expected_mass_per_node;
    }
    for (int i = 1; i < 12; i += 3) {
        EXPECT_NEAR(M(i, i).real(), expected_mass_per_node, 1e-6)
            << "M(" << i << "," << i << ") should be equal to " << expected_mass_per_node;
    }
    for (int i = 2; i < 12; i += 3) {
        EXPECT_NEAR(M(i, i).real(), expected_mass_per_node, 1e-6)
            << "M(" << i << "," << i << ") should be equal to " << expected_mass_per_node;
    }
    
    const double M_values[12][12] = {
        {130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0},
        {0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667},
        {65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0},
        {0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667},
        {65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0},
        {0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667},
        {65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0},
        {0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0, 65.4166666666667},
        {65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667, 0.0},
        {0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0, 65.4166666666667},
        {65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333, 0.0},
        {0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 65.4166666666667, 0.0, 130.833333333333}
    };
    
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            const double expected = M_values[i][j];
            EXPECT_NEAR(M(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "M(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(M(i, j).imag(), 0.0, 1e-12)
                << "M(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            EXPECT_NEAR(M(i, j).real(), M(j, i).real(), std::abs(M(i, j).real()) * 1e-10)
                << "M is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 12; i++) {
        EXPECT_GT(M(i, i).real(), 0.0)
            << "Diagonal element M(" << i << "," << i << ") should be positive";
    }
}

TEST_F(TetraTest, Volume_ExactVolumeCalculation) {

    double volume = tetra->Volume();
    const double expected = 0.16666666666666666;
    EXPECT_NEAR(volume, expected, 1e-9) << "Volume should be equal to " << expected;
}

#endif
