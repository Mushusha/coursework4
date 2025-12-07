// #define DISABLE_PYR_TESTS

#ifndef DISABLE_PYR_TESTS

#include "PyrTests.h"
#include <Eigen/Eigenvalues>

TEST_F(PyrTest, FF_ExactValues) {

    auto ff_node0 = pyr->FF(-1.0, -1.0, -1.0);
    EXPECT_NEAR(std::abs(ff_node0[0]), 1.0, 1e-9) << "FF[0] should be equal to 1.0";
    for (int i = 1; i < 5; i++) {
        EXPECT_NEAR(std::abs(ff_node0[i]), 0.0, 1e-9) << "FF[" << i << "] should be equal to 0.0";
    }

    auto ff_node1 = pyr->FF(1.0, -1.0, -1.0);
    EXPECT_NEAR(std::abs(ff_node1[1]), 1.0, 1e-9) << "FF[1] should be equal to 1.0";
    for (int i = 0; i < 5; i++) {
        if (i != 1) EXPECT_NEAR(std::abs(ff_node1[i]), 0.0, 1e-9) << "FF[" << i << "] should be equal to 0.0";
    }

    auto ff_apex = pyr->FF(0.0, 0.0, 1.0);
    EXPECT_NEAR(std::abs(ff_apex[4]), 1.0, 1e-9) << "FF[4] should be equal to 1.0";
    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(std::abs(ff_apex[i]), 0.0, 1e-9) << "FF[" << i << "] should be equal to 0.0";
    }

    auto ff_base_center = pyr->FF(0.0, 0.0, -1.0);
    const double expected_base = 0.25;
    const double expected_apex = 0.0;
    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(std::abs(ff_base_center[i]), expected_base, 1e-9)
            << "FF[" << i << "] should be equal to " << expected_base;
    }
    EXPECT_NEAR(std::abs(ff_base_center[4]), expected_apex, 1e-9)
        << "FF[4] should be equal to " << expected_apex;
}

TEST_F(PyrTest, B_ExactAllValues) {

    Eigen::MatrixXcd B = pyr->B(0.0, 0.0, 0.0);
    
    const double B_half = 0.5;
    const double B_values[6][15] = {
        {-B_half, B_half, B_half, -B_half, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {-B_half, -B_half, B_half, B_half, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, -B_half, B_half, B_half, -B_half, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -B_half, -B_half, B_half, B_half, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -B_half, B_half}
    };
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 15; j++) {
            const double expected = B_values[i][j];
            EXPECT_NEAR(B(i, j).real(), expected, std::abs(expected) * 1e-10)
                << "B(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(B(i, j).imag(), 0.0, 1e-12)
                << "B(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    EXPECT_EQ(B.rows(), 6) << "B should have 6 rows";
    EXPECT_EQ(B.cols(), 15) << "B should have 15 columns";
    
    double norm = B.norm();
    EXPECT_GT(norm, 0.0) << "Norm B should be positive";
    EXPECT_TRUE(std::isfinite(norm)) << "Norm B should be finite";
}

TEST_F(PyrTest, localK_ExactAllMatrixValues) {

    double E = 210e9;
    double nu = 0.3;
    pyr->set_constants(E, nu, 7850.0);
    
    double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    pyr->D = Eigen::MatrixXd::Zero(6, 6);
    pyr->D(0, 0) = pyr->D(1, 1) = pyr->D(2, 2) = factor * (1.0 - nu);
    pyr->D(0, 1) = pyr->D(0, 2) = pyr->D(1, 0) = pyr->D(1, 2) = pyr->D(2, 0) = pyr->D(2, 1) = factor * nu;
    pyr->D(3, 3) = pyr->D(4, 4) = pyr->D(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;
    
    Eigen::MatrixXcd K = pyr->localK();
    
    const double K_values[15][15] = {
        {9275521972.19254, 4111787192.84188, 2055893596.42094, -5163734779.35067, 822357438.568376, 411178719.284188, -5526911922.03823, -4111787192.84188, 411178719.284188, 3059839606.3331, -822357438.568376, 2055893596.42094, -1644714877.13675, 0, -4934144631.41026},
        {4111787192.84188, 9275521972.19254, 2055893596.42094, -822357438.568376, 3059839606.3331, 2055893596.42094, -4111787192.84188, -5526911922.03823, 411178719.284188, 822357438.568376, -5163734779.35067, 411178719.284188, 0, -1644714877.13675, -4934144631.41026},
        {2055893596.42094, 2055893596.42094, 5462032954.57918, -411178719.284188, 2055893596.42094, 705647834.683643, -411178719.284188, -411178719.284188, -1116826553.96783, 2055893596.42094, -411178719.284188, 705647834.683643, -3289429754.2735, -3289429754.2735, -5756502069.97863},
        {-5163734779.35067, -822357438.568376, -411178719.284188, 9275521972.19254, -4111787192.84188, -2055893596.42094, 3059839606.3331, 822357438.568376, -2055893596.42094, -5526911922.03823, 4111787192.84188, -411178719.284188, -1644714877.13675, 0, 4934144631.41026},
        {822357438.568376, 3059839606.3331, 2055893596.42094, -4111787192.84188, 9275521972.19254, 2055893596.42094, -822357438.568376, -5163734779.35067, 411178719.284188, 4111787192.84188, -5526911922.03823, 411178719.284188, 0, -1644714877.13675, -4934144631.41026},
        {411178719.284188, 2055893596.42094, 705647834.683643, -2055893596.42094, 2055893596.42094, 5462032954.57918, -2055893596.42094, -411178719.284188, 705647834.683642, 411178719.284188, -411178719.284188, -1116826553.96783, 3289429754.2735, -3289429754.2735, -5756502069.97863},
        {-5526911922.03823, -4111787192.84188, -411178719.284188, 3059839606.3331, -822357438.568376, -2055893596.42094, 9275521972.19254, 4111787192.84188, -2055893596.42094, -5163734779.35067, 822357438.568376, -411178719.284188, -1644714877.13675, 0, 4934144631.41026},
        {-4111787192.84188, -5526911922.03823, -411178719.284188, 822357438.568376, -5163734779.35067, -411178719.284188, 4111787192.84188, 9275521972.19254, -2055893596.42094, -822357438.568376, 3059839606.3331, -2055893596.42094, 0, -1644714877.13675, 4934144631.41026},
        {411178719.284188, 411178719.284188, -1116826553.96783, -2055893596.42094, 411178719.284188, 705647834.683642, -2055893596.42094, -2055893596.42094, 5462032954.57918, 411178719.284188, -2055893596.42094, 705647834.683642, 3289429754.2735, 3289429754.2735, -5756502069.97863},
        {3059839606.3331, 822357438.568376, 2055893596.42094, -5526911922.03823, 4111787192.84188, 411178719.284188, -5163734779.35067, -822357438.568376, 411178719.284188, 9275521972.19254, -4111787192.84188, 2055893596.42094, -1644714877.13675, 0, -4934144631.41026},
        {-822357438.568376, -5163734779.35067, -411178719.284188, 4111787192.84188, -5526911922.03823, -411178719.284188, 822357438.568376, 3059839606.3331, -2055893596.42094, -4111787192.84188, 9275521972.19254, -2055893596.42094, 0, -1644714877.13675, 4934144631.41026},
        {2055893596.42094, 411178719.284188, 705647834.683642, -411178719.284188, 411178719.284188, -1116826553.96783, -411178719.284188, -2055893596.42094, 705647834.683642, 2055893596.42094, -2055893596.42094, 5462032954.57918, -3289429754.2735, 3289429754.2735, -5756502069.97863},
        {-1644714877.13675, 0, -3289429754.2735, -1644714877.13675, 0, 3289429754.2735, -1644714877.13675, 0, 3289429754.2735, -1644714877.13675, 0, -3289429754.2735, 6578859508.54701, 0, 0},
        {0, -1644714877.13675, -3289429754.2735, 0, -1644714877.13675, -3289429754.2735, 0, -1644714877.13675, 3289429754.2735, 0, -1644714877.13675, 3289429754.2735, 0, 6578859508.54701, 0},
        {-4934144631.41026, -4934144631.41026, -5756502069.97863, 4934144631.41026, -4934144631.41026, -5756502069.97863, 4934144631.41026, 4934144631.41026, -5756502069.97863, -4934144631.41026, 4934144631.41026, -5756502069.97863, 0, 0, 23026008279.9145}
    };
    
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            const double expected = K_values[i][j];
            EXPECT_NEAR(K(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "K(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(K(i, j).imag(), 0.0, 1e-12)
                << "K(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            EXPECT_NEAR(K(i, j).real(), K(j, i).real(), std::abs(K(i, j).real()) * 1e-10)
                << "K is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 15; i++) {
        EXPECT_GT(K(i, i).real(), 0.0)
            << "Diagonal element K(" << i << "," << i << ") should be positive";
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(K.real());
    for (int i = 0; i < 15; i++) {
        EXPECT_GT(solver.eigenvalues()(i), 0.0)
            << "Eigenvalue " << i << " should be positive";
    }
}

TEST_F(PyrTest, localF_ExactAllValues) {

    double pressure_value = 1000.0;
    pyr->set_load(PRESSURE, 0, {pressure_value, 0, 0, 0, 0, 0});
    
    std::vector<double> F = pyr->localF();
    
    EXPECT_EQ(F.size(), 15) << "Vector size F should be 15";
    
    double total_force = 0.0;
    for (int i = 0; i < 15; i++) {
        total_force += std::abs(F[i]);
    }
    EXPECT_GT(total_force, 0.0) << "Total load should be positive";
    
    const double F_expected[15] = {
        0.0, 0.0, -250.0, 0.0, 0.0, -250.0, 0.0, 0.0, -250.0, 0.0, 0.0, -250.0, 0.0, 0.0, 0.0
    };
    
    for (int i = 0; i < 15; i++) {
        const double expected = F_expected[i];
        EXPECT_NEAR(F[i], expected, std::abs(expected) * 1e-6)
            << "F[" << i << "] should be equal to " << expected;
        EXPECT_TRUE(std::isfinite(F[i])) << "F[" << i << "] should be finite";
    }
}

TEST_F(PyrTest, localC_ExactAllMatrixValues) {

    Eigen::MatrixXd C = pyr->localC();
    
    const double C_values[5][5] = {
        {0.00235883034657922, 0.00193490277456276, 0.00172673158114712, 0.00193490277456276, 0.00399848090277778},
        {0.00193490277456276, 0.00235883034657922, 0.00193490277456276, 0.00172673158114712, 0.00399848090277778},
        {0.00172673158114712, 0.00193490277456276, 0.00235883034657922, 0.00193490277456276, 0.00399848090277778},
        {0.00193490277456276, 0.00172673158114712, 0.00193490277456276, 0.00235883034657922, 0.00399848090277778},
        {0.00399848090277778, 0.00399848090277778, 0.00399848090277778, 0.00399848090277778, 0.0176432291666667}
    };
    
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            const double expected = C_values[i][j];
            EXPECT_NEAR(C(i, j), expected, std::abs(expected) * 1e-10)
                << "C(" << i << "," << j << ") should be equal to " << expected;
        }
    }
    
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            EXPECT_NEAR(C(i, j), C(j, i), 1e-9)
                << "C is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            EXPECT_GT(C(i, j), 0.0)
                << "C(" << i << "," << j << ") should be positive";
        }
    }
}

TEST_F(PyrTest, localR_ExactAllValues) {

    std::vector<double> values(5, 1.0);
    std::vector<double> R = pyr->localR(values);
    
    const double expected = 0.06666666666666667;
    
    EXPECT_EQ(R.size(), 5) << "Vector size R should be 5";
    
    for (int i = 0; i < 5; i++) {
        EXPECT_NEAR(R[i], expected, std::abs(expected) * 1e-10)
            << "R[" << i << "] should be equal to " << expected;
        EXPECT_GT(R[i], 0.0) << "R[" << i << "] should be positive";
        EXPECT_TRUE(std::isfinite(R[i])) << "R[" << i << "] should be finite";
    }
    
    if (R.size() > 0) {
        double first_value = R[0];
        for (int i = 1; i < 5; i++) {
            EXPECT_NEAR(R[i], first_value, std::abs(first_value) * 1e-10)
                << "R[" << i << "] should be equal to R[0] = " << first_value;
        }
    }
}

TEST_F(PyrTest, localM_ExactAllMatrixValues) {

    Eigen::MatrixXcd M = pyr->localM();
    
    double volume = pyr->Volume();
    const double expected_mass_per_node = 7850.0 * volume / 5.0;
    
    for (int i = 0; i < 15; i += 3) {
        EXPECT_NEAR(M(i, i).real(), expected_mass_per_node, expected_mass_per_node * 0.2)
            << "M(" << i << "," << i << ") should be approximately " << expected_mass_per_node;
    }
    for (int i = 1; i < 15; i += 3) {
        EXPECT_NEAR(M(i, i).real(), expected_mass_per_node, expected_mass_per_node * 0.2)
            << "M(" << i << "," << i << ") should be approximately " << expected_mass_per_node;
    }
    for (int i = 2; i < 15; i += 3) {
        EXPECT_NEAR(M(i, i).real(), expected_mass_per_node, expected_mass_per_node * 0.2)
            << "M(" << i << "," << i << ") should be approximately " << expected_mass_per_node;
    }
    
    const double M_values[15][15] = {
        {18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 31.3880750868056, 0.0, 0.0},
        {0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 31.3880750868056, 0.0},
        {0.0, 0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 31.3880750868056},
        {15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 31.3880750868056, 0.0, 0.0},
        {0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 31.3880750868056, 0.0},
        {0.0, 0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 31.3880750868056},
        {13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 31.3880750868056, 0.0, 0.0},
        {0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 31.3880750868056, 0.0},
        {0.0, 0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 31.3880750868056},
        {15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 31.3880750868056, 0.0, 0.0},
        {0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 31.3880750868056, 0.0},
        {0.0, 0.0, 15.1889867803176, 0.0, 0.0, 13.5548429120049, 0.0, 0.0, 15.1889867803176, 0.0, 0.0, 18.5168182206469, 0.0, 0.0, 31.3880750868056},
        {31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 138.499348958333, 0.0, 0.0},
        {0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 138.499348958333, 0.0},
        {0.0, 0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 31.3880750868056, 0.0, 0.0, 138.499348958333}
    };
    
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            const double expected = M_values[i][j];
            EXPECT_NEAR(M(i, j).real(), expected, std::abs(expected) * 1e-6)
                << "M(" << i << "," << j << ") should be equal to " << expected;
            EXPECT_NEAR(M(i, j).imag(), 0.0, 1e-12)
                << "M(" << i << "," << j << ") should not have imaginary part";
        }
    }
    
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            EXPECT_NEAR(M(i, j).real(), M(j, i).real(), std::abs(M(i, j).real()) * 1e-10)
                << "M is not symmetric at position (" << i << "," << j << ")";
        }
    }
    
    for (int i = 0; i < 15; i++) {
        EXPECT_GT(M(i, i).real(), 0.0)
            << "Diagonal element M(" << i << "," << i << ") should be positive";
    }
}

TEST_F(PyrTest, Volume_ExactVolumeCalculation) {

    double vol = pyr->Volume();
    const double expected_volume = 0.3333333333333333;
    
    EXPECT_NEAR(vol, expected_volume, expected_volume * 0.1)
        << "Volume should be close to " << expected_volume << ", got " << vol;
    
    EXPECT_GT(vol, 0.0) << "Volume should be positive";
    EXPECT_TRUE(std::isfinite(vol)) << "Volume should be finite";
}

#endif
