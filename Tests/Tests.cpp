#include "Tests.h"

int Test() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}

// Node

TEST_F(NodeTest, ConstructorAndGetters) {
    EXPECT_EQ(node1->getID(), 1);
    EXPECT_DOUBLE_EQ(node1->getX(), 1.0);
    EXPECT_DOUBLE_EQ(node1->getY(), 2.0);
    EXPECT_DOUBLE_EQ(node1->getZ(), 3.0);

    EXPECT_EQ(node2->getID(), 2);
    EXPECT_DOUBLE_EQ(node2->getX(), 4.0);
    EXPECT_DOUBLE_EQ(node2->getY(), 5.0);
    EXPECT_DOUBLE_EQ(node2->getZ(), 6.0);
}

TEST_F(NodeTest, CopyConstructor) {
    Node nodeCopy(*node1);

    EXPECT_EQ(nodeCopy.getID(), 1);
    EXPECT_DOUBLE_EQ(nodeCopy.getX(), 1.0);
    EXPECT_DOUBLE_EQ(nodeCopy.getY(), 2.0);
    EXPECT_DOUBLE_EQ(nodeCopy.getZ(), 3.0);

    // Check that maps are copied
    node1->set_constraints(1, 10.5);
    node1->set_result(3.14, ResType::DISPLACEMENT);

    Node nodeCopy2(*node1);
    EXPECT_EQ(nodeCopy2.constraints.size(), 1);
    EXPECT_EQ(nodeCopy2.results.size(), 1);
}

TEST_F(NodeTest, CopyAssignment) {
    Node nodeCopy(0, { 0, 0, 0 });
    nodeCopy = *node1;

    EXPECT_EQ(nodeCopy.getID(), 1);
    EXPECT_DOUBLE_EQ(nodeCopy.getX(), 1.0);
    EXPECT_DOUBLE_EQ(nodeCopy.getY(), 2.0);
    EXPECT_DOUBLE_EQ(nodeCopy.getZ(), 3.0);

    nodeCopy = nodeCopy;
    EXPECT_EQ(nodeCopy.getID(), 1);
}

TEST_F(NodeTest, MoveConstructor) {
    Node nodeMoved(std::move(*node1));

    EXPECT_EQ(nodeMoved.getID(), 1);
    EXPECT_DOUBLE_EQ(nodeMoved.getX(), 1.0);
    EXPECT_DOUBLE_EQ(nodeMoved.getY(), 2.0);
    EXPECT_DOUBLE_EQ(nodeMoved.getZ(), 3.0);

    EXPECT_EQ(node1->getID(), 0);
    EXPECT_DOUBLE_EQ(node1->getX(), 0.0);
}

TEST_F(NodeTest, MoveAssignment) {
    Node nodeMoved(0, { 0, 0, 0 });
    nodeMoved = std::move(*node1);

    EXPECT_EQ(nodeMoved.getID(), 1);
    EXPECT_DOUBLE_EQ(nodeMoved.getX(), 1.0);
    EXPECT_DOUBLE_EQ(nodeMoved.getY(), 2.0);
    EXPECT_DOUBLE_EQ(nodeMoved.getZ(), 3.0);

    nodeMoved = std::move(nodeMoved);
    EXPECT_EQ(nodeMoved.getID(), 1);
}

TEST_F(NodeTest, SetConstraints) {
    node1->set_constraints(1, 10.5);
    node1->set_constraints(2, 20.5);

    EXPECT_EQ(node1->constraints.size(), 2);
    EXPECT_DOUBLE_EQ(node1->constraints[1], 10.5);
    EXPECT_DOUBLE_EQ(node1->constraints[2], 20.5);

    node1->set_constraints(1, 15.5);
    EXPECT_DOUBLE_EQ(node1->constraints[1], 15.5);
    EXPECT_EQ(node1->constraints.size(), 2);
}

TEST_F(NodeTest, SetResult) {
    node1->set_result(1.1, ResType::DISPLACEMENT);
    node1->set_result(2.2, ResType::STRESS);
    node1->set_result(3.3, ResType::DISPLACEMENT);

    EXPECT_EQ(node1->results.size(), 2);
    EXPECT_EQ(node1->results[ResType::DISPLACEMENT].size(), 2);
    EXPECT_EQ(node1->results[ResType::STRESS].size(), 1);

    EXPECT_DOUBLE_EQ(node1->results[ResType::DISPLACEMENT][0], 1.1);
    EXPECT_DOUBLE_EQ(node1->results[ResType::DISPLACEMENT][1], 3.3);
    EXPECT_DOUBLE_EQ(node1->results[ResType::STRESS][0], 2.2);
}

TEST_F(NodeTest, LoadMapOperations) {
    node1->load[1] = 100.0;
    node1->load[2] = 200.0;

    EXPECT_EQ(node1->load.size(), 2);
    EXPECT_DOUBLE_EQ(node1->load[1], 100.0);
    EXPECT_DOUBLE_EQ(node1->load[2], 200.0);

    Node nodeCopy(*node1);
    EXPECT_EQ(nodeCopy.load.size(), 2);

    Node nodeMoved(std::move(*node1));
    EXPECT_EQ(nodeMoved.load.size(), 2);
    EXPECT_EQ(node1->load.size(), 0);
}

// Elements
// Tri

TEST_F(TriTest, ConstructorAndBasicProperties) {
    EXPECT_EQ(tri->get_id(), 1);
    EXPECT_EQ(tri->get_type(), ElemType::TRI);
    EXPECT_EQ(tri->nodes_count(), 3);

    EXPECT_EQ(tri->get_node(0), 1);
    EXPECT_EQ(tri->get_node(1), 2);
    EXPECT_EQ(tri->get_node(2), 3);

    EXPECT_DOUBLE_EQ(tri->get_coord(0, 0), 0.0); // Node 1, x
    EXPECT_DOUBLE_EQ(tri->get_coord(1, 0), 1.0); // Node 2, x
    EXPECT_DOUBLE_EQ(tri->get_coord(2, 1), 1.0); // Node 3, y

    EXPECT_DOUBLE_EQ(tri->get_E(), 210e9);
    EXPECT_DOUBLE_EQ(tri->get_nu(), 0.3);
    EXPECT_DOUBLE_EQ(tri->get_rho(), 7850.0);
}

TEST_F(TriTest, GeometricProperties) {
    EXPECT_DOUBLE_EQ(tri->Volume(), expected_area);

    std::vector<double> inside_point = { 0.2, 0.2, 0.0 };
    std::vector<double> outside_point = { 1.0, 1.0, 0.0 };
    std::vector<double> vertex_point = { 0.0, 0.0, 0.0 };

    EXPECT_TRUE(tri->pointInElem(inside_point));
    EXPECT_FALSE(tri->pointInElem(outside_point));
    EXPECT_TRUE(tri->pointInElem(vertex_point));
}

TEST_F(TriTest, ShapeFunctions) {
    auto ff_node1 = tri->FF(0.3, 0.0, 0.0); 
    EXPECT_TRUE(ComplexNear(ff_node1[0], { 0.7, 0.0 }));
    EXPECT_TRUE(ComplexNear(ff_node1[1], { 0.3, 0.0 }));
    EXPECT_TRUE(ComplexNear(ff_node1[2], { 0.0, 0.0 }));

    auto ff_node2 = tri->FF(0.0, 0.3, 0.0); 
    EXPECT_TRUE(ComplexNear(ff_node2[0], { 0.7, 0.0 }));
    EXPECT_TRUE(ComplexNear(ff_node2[1], { 0.0, 0.0 }));
    EXPECT_TRUE(ComplexNear(ff_node2[2], { 0.3, 0.0 }));

    auto ff_center = tri->FF(1.0 / 3.0, 1.0 / 3.0, 0.0);
    for (const auto& val : ff_center) {
        EXPECT_NEAR(std::abs(val), 1.0 / 3.0, 1e-9);
    }
}

TEST_F(TriTest, BMatrix) {
    Eigen::MatrixXcd B = tri->B();

    EXPECT_EQ(B.rows(), 3);
    EXPECT_EQ(B.cols(), 6);

    EXPECT_NEAR(std::abs(B.coeff(0, 0).real()), 1.0, 1e-9);
    EXPECT_NEAR(std::abs(B.coeff(1, 1).real()), 1.0, 1e-9);
    EXPECT_NEAR(std::abs(B.coeff(2, 0).real()), 1.0, 1e-9);
}

TEST_F(TriTest, LocalStiffnessMatrix) {
    std::vector<double> x_coords = { 0.0, 1.0, 0.0 };
    std::vector<double> y_coords = { 0.0, 0.0, 1.0 };
    tri->set_coords(x_coords, y_coords, { 0.0, 0.0, 0.0 });

    double E = 210e9;
    double nu = 0.3;
    tri->set_constants(E, nu, 7850.0);

    double area = 0.5;
    double factor = E / (1 - nu * nu);

    Eigen::MatrixXd expected_D(3, 3);
    expected_D << 1, nu, 0,
        nu, 1, 0,
        0, 0, (1 - nu) / 2;
    expected_D *= factor;
    tri->D = expected_D;

    Eigen::MatrixXcd K = tri->localK();

    EXPECT_EQ(K.rows(), 6);
    EXPECT_EQ(K.cols(), 6);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(K(i, j).real(), K(j, i).real(), 1e-9);
            EXPECT_NEAR(K(i, j).imag(), 0.0, 1e-9); 
        }
    }

    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(K(i, i).real(), 0.0);
    }

    EXPECT_NEAR(K(0, 0).real(), K(1, 1).real(), 1e-9);
    EXPECT_NEAR(K(2, 2).real(), K(5, 5).real(), 1e-9);

    EXPECT_NEAR(K(0, 1).real() / 1e9, 75.0, 10.0);
}

TEST_F(TriTest, LocalMassMatrix) {
    std::vector<double> x_coords = { 0.0, 1.0, 0.0 };
    std::vector<double> y_coords = { 0.0, 0.0, 1.0 };
    tri->set_coords(x_coords, y_coords, { 0.0, 0.0, 0.0 });

    double density = 7850.0;
    tri->set_constants(210e9, 0.3, density);

    Eigen::MatrixXcd M = tri->localM();

    EXPECT_EQ(M.rows(), 6);
    EXPECT_EQ(M.cols(), 6);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(M(i, j).real(), M(j, i).real(), 1e-9);
            EXPECT_NEAR(M(i, j).imag(), 0.0, 1e-9);
        }
    }

    double expected_diag = density * expected_area / 6.0;
    double expected_off_diag = density * expected_area / 12.0;

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (i == j) {
                EXPECT_NEAR(M(i, j).real(), expected_diag, 1e-9);
            }
            else if ((i + j) % 2 == 0) { 
                EXPECT_NEAR(M(i, j).real(), expected_off_diag, 1e-9);
            }
            else {
                EXPECT_NEAR(M(i, j).real(), 0.0, 1e-9);
            }
        }
    }

    tri->set_constants(210e9, 0.3, 0.0);
    EXPECT_THROW(tri->localM(), std::runtime_error);
}

TEST_F(TriTest, LocalDampingMatrix) {
    std::vector<double> x_coords = { 0.0, 1.0, 0.0 };
    std::vector<double> y_coords = { 0.0, 0.0, 1.0 };
    tri->set_coords(x_coords, y_coords, { 0.0, 0.0, 0.0 });

    tri->set_constants(210e9, 0.3, 7850.0);

    Eigen::MatrixXd C = tri->localC();

    EXPECT_EQ(C.rows(), 3);
    EXPECT_EQ(C.cols(), 3);

    double expected_diag = expected_area / 6.0;
    double expected_off_diag = expected_area / 12.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) {
                EXPECT_NEAR(C(i, j), expected_diag, 1e-9);
            }
            else {
                EXPECT_NEAR(C(i, j), expected_off_diag, 1e-9);
            }
        }
    }

    EXPECT_GT(C.determinant(), 0.0);
}

TEST_F(TriTest, LoadAndForceVectors) {
    std::array<double, 6> pressure_value = { 1000.0, 0, 0, 0, 0, 0 };
    tri->set_load(PRESSURE, 0, pressure_value);

    std::vector<double> F = tri->localF();
    EXPECT_EQ(F.size(), 6);

    EXPECT_NEAR(F[1], 500.0, 1e-9);
    EXPECT_NEAR(F[3], 500.0, 1e-9); 

    std::vector<double> R = tri->localR({ 1.0, 1.0, 1.0 });
    EXPECT_EQ(R.size(), 3);
    EXPECT_DOUBLE_EQ(R[0], 1.0 * expected_area / 3.0);
}

TEST_F(TriTest, PressureLoadApplication) {
    std::array<double, 6> pressure_value = { 1000.0, 0, 0, 0, 0, 0 };
    tri->set_load(PRESSURE, 0, pressure_value);

    std::array<double, 6> pressure_value2 = { 2000.0, 0, 0, 0, 0, 0 };
    tri->set_load(PRESSURE, 1, pressure_value2);

    std::vector<double> F = tri->localF();

    EXPECT_NEAR(F[1], 500.0, 1e-9); 
    EXPECT_NEAR(F[3], -500.0, 1e-9);
    EXPECT_NEAR(F[4], -1000.0, 1e-9);
}

// Quad

TEST_F(QuadTest, ConstructorAndBasicProperties) {
    EXPECT_EQ(quad->get_id(), 1);
    EXPECT_EQ(quad->get_type(), ElemType::QUAD);
    EXPECT_EQ(quad->nodes_count(), 4);

    EXPECT_EQ(quad->get_node(0), 1);
    EXPECT_EQ(quad->get_node(1), 2);
    EXPECT_EQ(quad->get_node(3), 4);

    EXPECT_DOUBLE_EQ(quad->get_coord(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(quad->get_coord(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(quad->get_coord(2, 0), 1.0);
}

TEST_F(QuadTest, ShapeFunctions) {
    auto ff_node1 = quad->FF(-1.0, -1.0, 0.0);
    EXPECT_NEAR(std::abs(ff_node1[0]), 1.0, 1e-9);
    EXPECT_NEAR(std::abs(ff_node1[1]), 0.0, 1e-9);
    EXPECT_NEAR(std::abs(ff_node1[2]), 0.0, 1e-9);
    EXPECT_NEAR(std::abs(ff_node1[3]), 0.0, 1e-9);

    auto ff_center = quad->FF(0.0, 0.0, 0.0);
    for (const auto& val : ff_center) {
        EXPECT_NEAR(std::abs(val), 0.25, 1e-9);
    }
}

TEST_F(QuadTest, BMatrix) {
    Eigen::MatrixXcd B = quad->B(0.0, 0.0);

    EXPECT_EQ(B.rows(), 3);
    EXPECT_EQ(B.cols(), 8);

    EXPECT_NEAR(B(0, 0).real(), -0.5, 1e-9);
    EXPECT_NEAR(B(1, 1).real(), -0.5, 1e-9);
}

TEST_F(QuadTest, LocalStiffnessMatrix) {
    Eigen::MatrixXcd K = quad->localK();

    EXPECT_EQ(K.rows(), 8);
    EXPECT_EQ(K.cols(), 8);

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            EXPECT_NEAR(K(i, j).real(), K(j, i).real(), 1e-5);
            EXPECT_NEAR(K(i, j).imag(), 0.0, 1e-9);
        }
    }

    for (int i = 0; i < 8; ++i) {
        EXPECT_GT(K(i, i).real(), 0.0);
    }

    EXPECT_NEAR(K(0, 0).real() / 1e9, 103.0, 10.0);
}

TEST_F(QuadTest, LocalMassMatrix) {
    Eigen::MatrixXcd M = quad->localM();

    EXPECT_EQ(M.rows(), 8);
    EXPECT_EQ(M.cols(), 8);

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            EXPECT_NEAR(M(i, j).real(), M(j, i).real(), 1e-9);
            EXPECT_NEAR(M(i, j).imag(), 0.0, 1e-9);
        }
    }

    double expected_mass = 872.222;

    for (int i = 0; i < 6; i += 2) {
        EXPECT_NEAR(M(i, i).real(), expected_mass, expected_mass * 0.01);
        EXPECT_NEAR(M(i, i + 2).real(), expected_mass / 2.0, expected_mass * 0.01);
    }

    quad->set_constants(210e9, 0.3, 0.0);
    EXPECT_THROW(quad->localM(), std::runtime_error);
}

TEST_F(QuadTest, LocalDampingMatrix) {
    Eigen::MatrixXd C = quad->localC();

    EXPECT_EQ(C.rows(), 4);
    EXPECT_EQ(C.cols(), 4);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(C(i, j), C(j, i), 1e-9);
        }
    }

    for (int i = 0; i < 4; ++i) {
        EXPECT_GT(C(i, i), 0.0);
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j) {
                EXPECT_GT(C(i, j), C(i, (j + 1) % 4));
            }
        }
    }
}

TEST_F(QuadTest, LocalForceVector) {
    std::array<double, 6> pressure = { 1000.0, 0, 0, 0, 0, 0 };
    quad->set_load(PRESSURE, 0, pressure);

    std::vector<double> F = quad->localF();

    EXPECT_EQ(F.size(), 8);

    EXPECT_NEAR(F[1], 500.0, 1e-9); 
    EXPECT_NEAR(F[3], 500.0, 1e-9); 
    EXPECT_NEAR(F[5], 0.0, 1e-9);
    EXPECT_NEAR(F[7], 0.0, 1e-9);

    for (int i = 0; i < 8; i += 2) {
        EXPECT_NEAR(F[i], 0.0, 1e-9);
    }
}

TEST_F(QuadTest, LocalRVector) {
    std::vector<double> values = { 1.0, 2.0, 3.0, 4.0 };
    std::vector<double> R = quad->localR(values);

    EXPECT_EQ(R.size(), 4);

    for (int i = 0; i < 4; ++i) {
        EXPECT_GT(R[i], 0.0);
        EXPECT_NEAR(R[i], values[i] * 0.25, 0.01);
    }

    EXPECT_NEAR(R[1] / R[0], 2.0, 0.01);
    EXPECT_NEAR(R[2] / R[0], 3.0, 0.01);
}

TEST_F(QuadTest, PointInElement) {
    std::vector<double> inside_point = { 0.5, 0.5, 0.0 };
    std::vector<double> outside_point = { 1.5, 0.5, 0.0 };
    std::vector<double> vertex_point = { 0.0, 0.0, 0.0 };

    EXPECT_TRUE(quad->pointInElem(inside_point));
    EXPECT_FALSE(quad->pointInElem(outside_point));
    EXPECT_TRUE(quad->pointInElem(vertex_point));
}

TEST_F(QuadTest, VolumeCalculation) {
    EXPECT_NEAR(quad->Volume(), 1.0, 1e-9);
}

// infQuad

TEST_F(infQuadTest, StaticShapeFunctions) {
    inf_quad->is_dyn = false;

    auto ff = inf_quad->FF(-1.0 / 3.0, 0.0, 0.0);

    EXPECT_EQ(ff.size(), 4);

    EXPECT_NEAR(ff[0].real(), 0.375, 1e-9);  // Q0
    EXPECT_NEAR(ff[1].real(), 0.125, 1e-9); // C0
    EXPECT_NEAR(ff[2].real(), 0.125, 1e-9);  // C1
    EXPECT_NEAR(ff[3].real(), 0.375, 1e-9); // Q1

    for (const auto& val : ff) {
        EXPECT_NEAR(val.imag(), 0.0, 1e-9);
    }
}

TEST_F(infQuadTest, DynamicShapeFunctions) {
    inf_quad->is_dyn = true;

    auto ff = inf_quad->FF(-1.0 / 3.0, 0.0, 0.0);

    for (const auto& val : ff) {
        EXPECT_GT(std::abs(val.imag()), 0.0);
        EXPECT_TRUE(std::abs(val) > 0.0);
    }

    double original_omega = inf_quad->omega;
    inf_quad->omega = 200.0;
    auto ff_high_freq = inf_quad->FF(-1.0 / 3.0, 0.0, 0.0);
    EXPECT_NEAR(std::abs(ff[0]), std::abs(ff_high_freq[0]), 1e-9);

    inf_quad->omega = original_omega;
}