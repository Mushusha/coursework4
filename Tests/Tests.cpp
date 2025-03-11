#include "Tests.h"

double eps = 10e-06;

void tests() {
    test_triElement_B();
    test_triElement_locK();
    test_triElement_globK();
    test_triElement_locF();
    test_triElement_globF();
    test_triElement_constraintsK();
    test_triElement_constraintsF();
    test_triElement_dispToElem();
    test_triElement_calcStrain();
    test_triElement_calcStress();
}

// TRIANGLE TESTS

void test_triElement_B() {
    triElement elem;
    std::vector<double> x = { 1, 4, 9 };
    std::vector<double> y = { 1, 6, 1 };
    std::vector<double> z = { 0, 0, 0 };

    elem.set_coords(x, y, z);

    Eigen::MatrixXd expectedB(3, 6);
    expectedB <<
        0.125, 0, 0, 0, -0.125, 0,
        0, 0.125, 0, -0.2, 0, 0.075,
        0.125, 0.125, -0.2, 0, 0.075, -0.125;
    Eigen::MatrixXd resultB = elem.B();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 6; ++j)
            assert(resultB(i, j) == expectedB(i, j));

    std::cout << "triElement_B: Test passed!" << std::endl;
}

void test_triElement_locK() {
    triElement elem;
    std::vector<double> x = { 1, 4, 9 };
    std::vector<double> y = { 1, 6, 1 };
    std::vector<double> z = { 0, 0, 0 };

    elem.set_coords(x, y, z);
    elem.set_constants(1, 0.6, 1);

    Eigen::MatrixXd expectedlocK(6, 6);
    expectedlocK <<
        -0.29296875, -0.48828125, -0.15625, 0.9375, 0.44921875, -0.44921875,
        -0.48828125, -0.29296875, -0.15625, 0.625, 0.64453125, -0.33203125,
        -0.15625, -0.15625, 0.25, 0, -0.09375, 0.15625,
        0.9375, 0.625, 0, -1, -0.9375, 0.375,
        0.44921875, 0.64453125, -0.09375, -0.9375, -0.35546875, 0.29296875,
        -0.44921875, -0.33203125, 0.15625, 0.375, 0.29296875, -0.04296875;

    Eigen::MatrixXd resultlocK = elem.localK();
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            assert(std::abs(resultlocK(i, j) - expectedlocK(i, j)) < eps);

    std::cout << "triElement_locK: Test passed!" << std::endl;
}

void test_triElement_globK() {
    string file = "test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalK();
    Eigen::SparseMatrix<double> resultGlobK = data.K;

    Eigen::MatrixXd expectedGlobK = Eigen::MatrixXd::Zero(10, 10);
    expectedGlobK <<
        1.209065363, 0.8614259729, -0.490018733, 0.2678571429, 0, 0, 0.4437730467, -0.2364259729, -1.162819677, -0.8928571429,
        0.8614259729, 1.208919808, -0.2678571429, 0.4049057915, 0, 0, 0.2992883128, -0.4468711589, -0.8928571429, -1.16695444,
        -0.490018733, -0.2678571429, 1.209065363, -0.8614259729, 0.4437730467, 0.2364259729, 0, 0, -1.162819677, 0.8928571429,
        0.2678571429, 0.4049057915, -0.8614259729, 1.208919808, -0.2992883128, -0.4468711589, 0, 0, 0.8928571429, -1.16695444,
        0, 0, 0.4437730467, -0.2992883128, 1.297135271, 0.9242883128, -0.4019488253, 0.2678571429, -1.338959492, -0.8928571429,
        0, 0, 0.2364259729, -0.4468711589, 0.9242883128, 1.297302875, -0.2678571429, 0.4932888588, -0.8928571429, -1.343720575,
        0.4437730467, 0.2992883128, 0, 0, -0.4019488253, -0.2678571429, 1.297135271, -0.9242883128, -1.338959492, 0.8928571429,
        -0.2364259729, -0.4468711589, 0, 0, 0.2678571429, 0.4932888588, -0.9242883128, 1.297302875, 0.8928571429, -1.343720575,
        -1.162819677, -0.8928571429, -1.162819677, 0.8928571429, -1.338959492, -0.8928571429, -1.338959492, 0.8928571429, 5.003558338, 0,
        -0.8928571429, -1.16695444, 0.8928571429, -1.16695444, -0.8928571429, -1.343720575, 0.8928571429, -1.343720575, 0, 5.02135003;
    //for (int elem = 0; elem < 4; elem++) {
    //    Eigen::MatrixXd k(6, 6);
    //    k = data.getElem(elem)->localK();
    //    for (int i = 0; i < 3; i++)
    //        for (int j = 0; j < 3; j++) {
    //            expectedGlobK(2 * data.getElem(elem)->get_nodes(i) - 2, 2 * data.getElem(elem)->get_nodes(j) - 2) += k(2 * i, 2 * j);
    //            expectedGlobK(2 * data.getElem(elem)->get_nodes(i) - 1, 2 * data.getElem(elem)->get_nodes(j) - 2) += k(2 * i + 1, 2 * j);
    //            expectedGlobK(2 * data.getElem(elem)->get_nodes(i) - 2, 2 * data.getElem(elem)->get_nodes(j) - 1) += k(2 * i, 2 * j + 1);
    //            expectedGlobK(2 * data.getElem(elem)->get_nodes(i) - 1, 2 * data.getElem(elem)->get_nodes(j) - 1) += k(2 * i + 1, 2 * j + 1);
    //        }
    //}

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; ++j)
            assert(std::abs(resultGlobK.coeffRef(i, j) - expectedGlobK(i, j)) < eps);

    std::cout << "triElement_GlobK: Test passed!" << std::endl;
}

void test_triElement_locF() {
    triElement elem;
    std::vector<double> x = { 1, 1, 5 };
    std::vector<double> y = { 1, 4, 1 };
    std::vector<double> z = { 0, 0, 0 };

    elem.set_coords(x, y, z);
    elem.set_load(4, 0, std::array<double, 6>{2, 0, 0, 0, 0, 0});
    elem.set_load(4, 2, std::array<double, 6>{3, 0, 0, 0, 0 ,0});

    std::vector<double> expectedlocF{3, 6, 3, 0, 0, 6};

    std::vector<double> resultlocF = elem.localF();

    for (int i = 0; i < 6; i++)
        assert(std::abs(resultlocF[i] - expectedlocF[i]) < eps);

    std::cout << "triElement_locF: Test passed!" << std::endl;
}

void test_triElement_globF() {
    string file = "test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalF();
    Eigen::SparseVector<double> resultGlobF = data.F;

    Eigen::SparseVector<double> expectedGlobF;
    expectedGlobF.resize(10);

    expectedGlobF.insert(1) = -1;
    expectedGlobF.insert(3) = -1;

    for (int i = 0; i < 10; i++)
        assert(std::abs(resultGlobF.coeffRef(i) - expectedGlobF.coeffRef(i)) < eps);

    std::cout << "triElement_GlobF: Test passed!" << std::endl;

    assert("xyu" == "zalypa");
}


void test_triElement_constraintsK() {
    string file = "test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalK();
    data.fillGlobalF();
    data.fillConstraints();
    Eigen::SparseMatrix<double> resultGlobK = data.K;

    Eigen::MatrixXd expectedGlobK = Eigen::MatrixXd::Zero(10, 10);
    expectedGlobK <<
        1.209065363, 0.8614259729, -0.490018733, 0.2678571429, 0, 0, 0.4437730467, 0, -1.162819677, -0.8928571429,
        0.8614259729, 1.208919808, -0.2678571429, 0.4049057915, 0, 0, 0.2992883128, 0, -0.8928571429, -1.16695444,
        -0.490018733, -0.2678571429, 1.209065363, -0.8614259729, 0.4437730467, 0, 0, 0, -1.162819677, 0.8928571429,
        0.2678571429, 0.4049057915, -0.8614259729, 1.208919808, -0.2992883128, 0, 0, 0, 0.8928571429, -1.16695444,
        0, 0, 0.4437730467, -0.2992883128, 1.297135271, 0, -0.4019488253, 0, -1.338959492, -0.8928571429,
        0, 0, 0, 0, 0, 1.297302875, 0, 0, 0, 0,
        0.4437730467, 0.2992883128, 0, 0, -0.4019488253, 0, 1.297135271, 0, -1.338959492, 0.8928571429,
        0, 0, 0, 0, 0, 0, 0, 1.297302875, 0, 0,
        -1.162819677, -0.8928571429, -1.162819677, 0.8928571429, -1.338959492, 0, -1.338959492, 0, 5.003558338, 0,
        -0.8928571429, -1.16695444, 0.8928571429, -1.16695444, -0.8928571429, 0, 0.8928571429, 0, 0, 5.02135003;

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            assert(std::abs(resultGlobK.coeffRef(i, j) - expectedGlobK(i, j)) < eps);

    std::cout << "triElement_constraintsK: Test passed!" << std::endl;
}

void test_triElement_constraintsF() {
    string file = "test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalK();
    data.fillGlobalF();
    data.fillConstraints();
    Eigen::SparseVector<double> resultGlobF = data.F;

    Eigen::SparseVector<double> expectedGlobF;
    expectedGlobF.resize(10);

    expectedGlobF.insert(1) = -1;
    expectedGlobF.insert(3) = -1;

    for (int i = 0; i < 10; i++)
        assert(std::abs(resultGlobF.coeffRef(i) - expectedGlobF.coeffRef(i)) < eps);

    std::cout << "triElement_constraintsR: Test passed!" << std::endl;
}

void test_triElement_dispToElem() {
    string file = "test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.solve();
    data.U.coeffRef(0) = 0.5600011511; // 1
    data.U.coeffRef(1) = -1.679999998;
    data.U.coeffRef(2) = -0.5599988481; // 2
    data.U.coeffRef(3) = -1.679999998;
    data.U.coeffRef(4) = -0.5599988478; // 3
    data.U.coeffRef(5) = 0;
    data.U.coeffRef(6) = 0.5600011508; // 4
    data.U.coeffRef(7) = 0;
    data.U.coeffRef(8) = 1.151483577e-06; // 5
    data.U.coeffRef(9) = -0.7808591097;

    //data.fillFields();

    std::array <std::vector <double>, 4> resultU;
    //for (int i = 0; i < 4; i++)
    //    resultU[i] = data.get_elem(i)->results[DISPLACEMENT];

    std::array <std::vector <double>, 4> expectedU;
    expectedU[0] = std::vector<double>{ -0.5599988481, -1.679999998, 1.151483577e-06, -0.7808591097, 0.5600011511, -1.679999998 };
    expectedU[1] = std::vector<double>{ 0.5600011508, 0, 1.151483577e-06, -0.7808591097, -0.5599988478, 0 };
    expectedU[2] = std::vector<double>{ 1.151483577e-06, -0.7808591097, 0.5600011508, 0, 0.5600011511, -1.679999998 };
    expectedU[3] = std::vector<double>{ 1.151483577e-06, -0.7808591097, -0.5599988481, -1.679999998, -0.5599988478, 0 };

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 6; j++)
            assert(std::abs(resultU[i][j] - expectedU[i][j]) < eps);

     std::cout << "triElement_dispToElem: Test passed!" << std::endl;
}

void test_triElement_calcStrain() {
    string file = "test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.solve();
    data.U.coeffRef(0) = 0.5600011511; // 1
    data.U.coeffRef(1) = -1.679999998;
    data.U.coeffRef(2) = -0.5599988481; // 2
    data.U.coeffRef(3) = -1.679999998;
    data.U.coeffRef(4) = -0.5599988478; // 3
    data.U.coeffRef(5) = 0;
    data.U.coeffRef(6) = 0.5600011508; // 4
    data.U.coeffRef(7) = 0;
    data.U.coeffRef(8) = 1.151483577e-06; // 5
    data.U.coeffRef(9) = -0.7808591097;

    //data.fillFields();

    std::array <std::vector <double>, 4> resultEpsilon;
    //for (int i = 0; i < 4; i++)
    //    resultEpsilon[i] = data.get_elem(i)->results[STRAIN];

    std::array <std::vector <double>, 4> expectedEpsilon;
    expectedEpsilon[0] = std::vector<double>{ 0.5599999996, -0.8400015764, 0.0 };
    expectedEpsilon[1] = std::vector<double>{ 0.5599999993, -0.8400013787, 0.0 };
    expectedEpsilon[2] = std::vector<double>{ 0.5599999995, -0.8399999999, 0.0000001508 };
    expectedEpsilon[3] = std::vector<double>{ 0.5599999994, -0.8399999999, 0.0000001508 };

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
            assert(std::abs(resultEpsilon[i][j] - expectedEpsilon[i][j]) < eps);

    std::cout << "triElement_calcStrain: Test passed!" << std::endl;
}

void test_triElement_calcStress() {
    string file = "test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.solve();
    data.U.coeffRef(0) = 0.5600011511; // 1
    data.U.coeffRef(1) = -1.679999998;
    data.U.coeffRef(2) = -0.5599988481; // 2
    data.U.coeffRef(3) = -1.679999998;
    data.U.coeffRef(4) = -0.5599988478; // 3
    data.U.coeffRef(5) = 0;
    data.U.coeffRef(6) = 0.5600011508; // 4
    data.U.coeffRef(7) = 0;
    data.U.coeffRef(8) = 1.151483577e-06; // 5
    data.U.coeffRef(9) = -0.7808591097;

    //data.fillFields();

    std::array <std::vector <double>, 4> resultSigma;
    //for (int i = 0; i < 4; i++)
    //    resultSigma[i] = data.get_elem(i)->results[STRESS];

    std::array <std::vector <double>, 4> expectedSigma;
    expectedSigma[0] = std::vector<double>{ 0.0000005471, -1.0000065786, 0 };
    expectedSigma[1] = std::vector<double>{ 0.0000008289, -1.0000061554, 0 };
    expectedSigma[2] = std::vector<double>{ 0.0000027991, -1.0000032005, 0.0000000539 };
    expectedSigma[3] = std::vector<double>{ 0.0000027991, -1.0000032006, 0.0000000539 };

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
            assert(std::abs(resultSigma[i][j] - expectedSigma[i][j]) < eps);

    std::cout << "triElement_calcStress: Test passed!" << std::endl;
}
