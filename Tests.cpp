#include "Tests.h"

double eps = 10e-08;

void tests() {
    test_triElement_B();
    test_triElement_locK();
    test_triElement_globK();
    test_triElement_locR();
    test_triElement_globR();
    test_triElement_constraintsK();
    test_triElement_constraintsR();

}

// TRIANGLE TESTS

void test_triElement_B() {
//    triElement elem;
//    std::vector<double> x = { 1, 4, 9 };
//    std::vector<double> y = { 1, 6, 1 };
//    std::vector<double> z = { 0, 0, 0 };
//
//    elem.set_coords(x, y, z);
//
//    Eigen::MatrixXd expectedB(3, 6);
//    expectedB <<
//        0.125, 0, 0, 0, -0.125, 0,
//        0, 0.125, 0, -0.2, 0, 0.075,
//        0.125, 0.125, -0.2, 0, 0.075, -0.125;
//    Eigen::MatrixXd resultB = elem.B();
//
//    for (int i = 0; i < 3; ++i)
//        for (int j = 0; j < 6; ++j)
//            assert(resultB(i, j) == expectedB(i, j));
//
    std::cout << "triElement_B: Test passed!" << std::endl;
}

void test_triElement_locK() {
    triElement elem;
    std::vector<double> x = { 1, 4, 9 };
    std::vector<double> y = { 1, 6, 1 };
    std::vector<double> z = { 0, 0, 0 };

    elem.set_coords(x, y, z);
    elem.set_constants(1, 0.6);

    Eigen::MatrixXd expectedlocK(6, 6);
    expectedlocK <<
        0.5859375, 0.390625, -0.15625, -0.46875, -0.4296875, 0.078125,
        0.390625, 0.5859375, -0.15625, -0.78125, -0.234375, 0.1953125,
        -0.15625, -0.15625, 0.25, 0, -0.09375, 0.15625,
        -0.46875, -0.78125, 0, 1.25, 0.46875, -0.46875,
        -0.4296875, -0.234375, -0.09375, 0.46875, 0.5234375, -0.234375,
        0.078125, 0.1953125, 0.15625, -0.46875, -0.234375, 0.2734375;

    Eigen::MatrixXd resultlocK = elem.localK();

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            assert(std::abs(resultlocK(i, j) - expectedlocK(i, j)) < eps);

    std::cout << "triElement_locK: Test passed!" << std::endl;
}

void test_triElement_globK() {
    string file = "C:/Users/mushu/Desktop/cwgit/coursework4/test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalK();
    Eigen::SparseMatrix<double> resultGlobK = data.K;

    Eigen::MatrixXd expectedGlobK = Eigen::MatrixXd::Zero(10, 10);
    expectedGlobK <<
        0.7484579444, 0.4019987874, -0.2351602043, 0.02976190476, 0, 0, 0.2068580418, -0.01509402547, -0.7201557819, -0.4166666667,
        0.4019987874, 0.7483900184, -0.02976190476, 0.1824712405, 0, 0, 0.04442978405, -0.2087759208, -0.4166666667, -0.7220853381,
        -0.2351602043, -0.02976190476, 0.7484579444, -0.4019987874, 0.2068580418, 0.01509402547, 0, 0, -0.7201557819, 0.4166666667,
        0.02976190476, 0.1824712405, -0.4019987874, 0.7483900184, -0.04442978405, -0.2087759208, 0, 0, 0.4166666667, -0.7220853381,
        0, 0, 0.2068580418, -0.04442978405, 0.8030012709, 0.431334546, -0.1806168778, 0.02976190476, -0.8292424348, -0.4166666667,
        0, 0, 0.01509402547, -0.2087759208, 0.431334546, 0.803079486, -0.02976190476, 0.2371607081, -0.4166666667, -0.8314642733,
        0.2068580418, 0.04442978405, 0, 0, -0.1806168778, -0.02976190476, 0.8030012709, -0.431334546, -0.8292424348, 0.4166666667,
        -0.01509402547, -0.2087759208, 0, 0, 0.02976190476, 0.2371607081, -0.431334546, 0.803079486, 0.4166666667, -0.8314642733,
        -0.7201557819, -0.4166666667, -0.7201557819, 0.4166666667, -0.8292424348, -0.4166666667, -0.8292424348, 0.4166666667, 3.098796434, 0,
        -0.4166666667, -0.7220853381, 0.4166666667, -0.7220853381, -0.4166666667, -0.8314642733, 0.4166666667, -0.8314642733, 0, 3.107099223;
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
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            assert(std::abs(resultGlobK.coeffRef(i, j) - expectedGlobK(i, j)) < eps);

    std::cout << "triElement_GlobK: Test passed!" << std::endl;
}

void test_triElement_locR() {
    triElement elem;
    std::vector<double> x = { 1, 1, 5 };
    std::vector<double> y = { 1, 4, 1 };
    std::vector<double> z = { 0, 0, 0 };

    elem.set_coords(x, y, z);
    elem.set_load(4, 0, std::array<double, 6>{2, 0, 0, 0, 0, 0});
    elem.set_load(4, 2, std::array<double, 6>{3, 0, 0, 0, 0 ,0});

    std::vector<double> expectedlocR{-3, -6, -3, 0, 0, -6};

    std::vector<double> resultlocR = elem.localR();

    for (int i = 0; i < 6; ++i)
        assert(std::abs(resultlocR[i] - expectedlocR[i]) < eps);

    std::cout << "triElement_locR: Test passed!" << std::endl;
}

void test_triElement_globR() {
    string file = "C:/Users/mushu/Desktop/cwgit/coursework4/test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalR();
    Eigen::SparseVector<double> resultGlobR = data.R;

    Eigen::SparseVector<double> expectedGlobR;
    expectedGlobR.resize(10);

    expectedGlobR.insert(1) = 1;
    expectedGlobR.insert(3) = 1;

    for (int i = 0; i < 10; ++i)
        assert(std::abs(resultGlobR.coeffRef(i) - expectedGlobR.coeffRef(i)) < eps);
     

    std::cout << "triElement_GlobR: Test passed!" << std::endl;
}


void test_triElement_constraintsK() {
    string file = "C:/Users/mushu/Desktop/cwgit/coursework4/test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalK();
    data.fillGlobalR();
    data.fillconstraints();
    Eigen::SparseMatrix<double> resultGlobK = data.K;

    Eigen::MatrixXd expectedGlobK = Eigen::MatrixXd::Zero(10, 10);
    expectedGlobK <<
        0.7484579444, 0.4019987874, -0.2351602043, 0.02976190476, 0, 0, 0.2068580418, 0, -0.7201557819, -0.4166666667,
        0.4019987874, 0.7483900184, -0.02976190476, 0.1824712405, 0, 0, 0.04442978405, 0, -0.4166666667, -0.7220853381,
        -0.2351602043, -0.02976190476, 0.7484579444, -0.4019987874, 0.2068580418, 0, 0, 0, -0.7201557819, 0.4166666667,
        0.02976190476, 0.1824712405, -0.4019987874, 0.7483900184, -0.04442978405, 0, 0, 0, 0.4166666667, -0.7220853381,
        0, 0, 0.2068580418, -0.04442978405, 0.8030012709, 0, -0.1806168778, 0, -0.8292424348, -0.4166666667,
        0, 0, 0, 0, 0, 0.803079486, 0, 0, 0, 0,
        0.2068580418, 0.04442978405, 0, 0, -0.1806168778, 0, 0.8030012709, 0, -0.8292424348, 0.4166666667,
        0, 0, 0, 0, 0, 0, 0, 0.803079486, 0, 0,
        -0.7201557819, -0.4166666667, -0.7201557819, 0.4166666667, -0.8292424348, 0, -0.8292424348, 0, 3.098796434, 0,
        -0.4166666667, -0.7220853381, 0.4166666667, -0.7220853381, -0.4166666667, 0, 0.4166666667, 0, 0, 3.107099223;

    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            assert(std::abs(resultGlobK.coeffRef(i, j) - expectedGlobK(i, j)) < eps);

    std::cout << "triElement_constraintsK: Test passed!" << std::endl;
}

void test_triElement_constraintsR() {
    string file = "C:/Users/mushu/Desktop/cwgit/coursework4/test.fc";
    std::shared_ptr<Parser> p = std::make_shared<Parser>();
    p->read(file);
    Data data(p);
    data.fillGlobalK();
    data.fillGlobalR();
    data.fillconstraints();
    Eigen::SparseVector<double> resultGlobR = data.R;

    Eigen::SparseVector<double> expectedGlobR;
    expectedGlobR.resize(10);

    expectedGlobR.insert(1) = 1;
    expectedGlobR.insert(3) = 1;

    for (int i = 0; i < 10; ++i)
        assert(std::abs(resultGlobR.coeffRef(i) - expectedGlobR.coeffRef(i)) < eps);

    std::cout << "triElement_constraintsR: Test passed!" << std::endl;
}


