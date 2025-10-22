#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cassert>

#include "Eigen/Core"
#include "Element.h"
#include "MathMV.h"

template <int NODES>
class SpectralQuad : public Element {
    static_assert(NODES >= 2, "NODES must be >= 2");

public:
    SpectralQuad() : Element() {
        init_gll(); 
    }
    SpectralQuad(int id, ElemType type, std::vector<int> nodes)
        : Element(id, type, nodes) {
        init_gll();
    }
    SpectralQuad(const SpectralQuad& other) : Element(other) {
        init_gll(); 
    }
    SpectralQuad& operator=(const SpectralQuad& other) {
        if (this != &other) Element::operator=(other);
        return *this;
    }
    SpectralQuad(SpectralQuad&& other) noexcept : Element(std::move(other)) {
        init_gll(); 
    }
    SpectralQuad& operator=(SpectralQuad&& other) noexcept {
        if (this != &other) Element::operator=(std::move(other));
        return *this;
    }
    virtual ~SpectralQuad() = default;


    std::vector<std::complex<double>> FF(double ksi, double eta, double zeta = 0) override;

    Eigen::MatrixXcd B(double ksi = 0, double eta = 0, double zeta = 0) override;

    Eigen::MatrixXcd localK() override;
    std::vector<double> localF(double mult = 1) override;

    Eigen::MatrixXd localC() override;
    std::vector<double> localR(std::vector<double> value) override;
    Eigen::MatrixXcd localM() override;

    std::vector<int> edge_to_node(int edge) final;

    double gaussPoint(LocVar var, int i) override;
    double weight(LocVar var, int i) override;
    double Volume() final;

    Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta = 0) override;
    Eigen::MatrixXcd J(double ksi, double eta, double zeta = 0) override;

    bool pointInElem(std::vector<double> point) override;
    std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;

protected:
    void set_pressure(int edge, double value) override;

private:
    static constexpr int nNodes = NODES * NODES;
    std::vector<double> gll_x;
    std::vector<double> gll_w;

    void init_gll();

    double lagrange1D(int i, double x) const;
    double dlagrange1D(int i, double x) const;

    double len_edge(int edge);
};


template<int NODES>
void SpectralQuad<NODES>::init_gll() {
    if (gll_x.empty() || gll_w.empty()) 
        compute_gll_nodes_weights(NODES - 1, gll_x, gll_w);
}

template<int NODES>
double SpectralQuad<NODES>::lagrange1D(int i, double x) const {
    double val = 1.0;
    for (int m = 0; m < NODES; ++m) {
        if (m == i) continue;
        val *= (x - gll_x[m]) / (gll_x[i] - gll_x[m]);
    }
    return val;
}

template<int NODES>
double SpectralQuad<NODES>::dlagrange1D(int i, double x) const {
    double sum = 0.0;
    for (int m = 0; m < NODES; ++m) {
        if (m == i) 
            continue;

        double prod = 1.0 / (gll_x[i] - gll_x[m]);

        for (int k = 0; k < NODES; ++k) {
            if (k == i || k == m) 
                continue;

            prod *= (x - gll_x[k]) / (gll_x[i] - gll_x[k]);
        }
        sum += prod;
    }
    return sum;
}

template<int NODES>
std::vector<std::complex<double>> SpectralQuad<NODES>::FF(double ksi, double eta, double) {
    std::vector<std::complex<double>> N(nNodes);
    int idx = 0;

    for (int j = 0; j < NODES; ++j) {
        double Nj_eta = lagrange1D(j, eta);
        for (int i = 0; i < NODES; ++i) {
            double Ni_ksi = lagrange1D(i, ksi);
            double value = Ni_ksi * Nj_eta;
            N[idx] = std::complex<double>(value, 0.0);

            idx++;
        }
    }
    return N;
}

template<int NODES>
Eigen::MatrixXcd SpectralQuad<NODES>::gradFF(double ksi, double eta, double) {
    Eigen::MatrixXcd grad(2, nNodes);
    grad.setZero();
    int idx = 0;

    for (int j = 0; j < NODES; ++j) {
        for (int i = 0; i < NODES; ++i, ++idx) {
            double dN_dksi = dlagrange1D(i, ksi) * lagrange1D(j, eta);
            double dN_deta = lagrange1D(i, ksi) * dlagrange1D(j, eta);

            grad(0, idx) = std::complex<double>(dN_dksi, 0.0);
            grad(1, idx) = std::complex<double>(dN_deta, 0.0);
        }
    }
    return grad;
}

template<int NODES>
Eigen::MatrixXcd SpectralQuad<NODES>::J(double ksi, double eta, double) {
    Eigen::MatrixXcd jac = Eigen::MatrixXcd::Zero(2, 2);
    Eigen::MatrixXcd g = gradFF(ksi, eta);

    for (int a = 0; a < nNodes; ++a) {
        jac(0, 0) += g(0, a) * x[a];
        jac(0, 1) += g(0, a) * y[a];
        jac(1, 0) += g(1, a) * x[a];
        jac(1, 1) += g(1, a) * y[a];
    }
    return jac;
}

template<int NODES>
Eigen::MatrixXcd SpectralQuad<NODES>::B(double ksi, double eta, double) {
    Eigen::MatrixXcd Bm = Eigen::MatrixXcd::Zero(3, 2 * nNodes);
    Eigen::Matrix2cd invJ = J(ksi, eta).inverse();
    Eigen::MatrixXcd dN = invJ * gradFF(ksi, eta);

    for (int a = 0; a < nNodes; ++a) {
        double dNx = dN(0, a).real();
        double dNy = dN(1, a).real();

        Bm(0, 2 * a) = dNx;
        Bm(1, 2 * a + 1) = dNy;
        Bm(2, 2 * a) = dNy;
        Bm(2, 2 * a + 1) = dNx;
    }
    return Bm;
}

template<int NODES>
Eigen::MatrixXcd SpectralQuad<NODES>::localK() {
    Eigen::MatrixXcd K = Eigen::MatrixXcd::Zero(2 * nNodes, 2 * nNodes);

    for (int j = 0; j < NODES; ++j) {
        double eta = gll_x[j];
        double w_eta = gll_w[j];

        for (int i = 0; i < NODES; ++i) {
            double ksi = gll_x[i];
            double w_ksi = gll_w[i];

            Eigen::MatrixXcd Bm = B(ksi, eta);
            Eigen::MatrixXcd Jm = J(ksi, eta);

            double detJ = std::abs(Jm.determinant());

            K += (w_ksi * w_eta) * (Bm.transpose() * D * Bm * detJ);
        }
    }
    return K;
}

template<int NODES>
std::vector<double> SpectralQuad<NODES>::localF(double mult) {
    std::vector<double> F(2 * nNodes, 0.0);
    if (load.empty()) return F;

    for (auto const& l : load) {
        int edge = l.first.first;
        int comp = l.first.second;
        double value = l.second;

        std::vector<int> edge_nodes = edge_to_node(edge);
        int m = edge_nodes.size();

        for (int i = 0; i < m; ++i) {
            double ksi, eta;
            double weight = gll_w[i];

            if (edge == 0) {
                ksi = gll_x[i];
                eta = -1.0;
            }
            else if (edge == 1) {
                ksi = 1.0;
                eta = gll_x[i];
            }
            else if (edge == 2) {
                ksi = gll_x[i];
                eta = 1.0;
            }
            else {
                ksi = -1.0;
                eta = gll_x[i];
            }

            Eigen::MatrixXcd grad = gradFF(ksi, eta);
            double dx = 0.0, dy = 0.0;

            if (edge == 0 || edge == 2) {
                for (int a = 0; a < nNodes; ++a) {
                    dx += grad(0, a).real() * x[a];
                    dy += grad(0, a).real() * y[a];
                }
            }
            else {
                for (int a = 0; a < nNodes; ++a) {
                    dx += grad(1, a).real() * x[a];
                    dy += grad(1, a).real() * y[a];
                }
            }

            double detJ_edge = std::sqrt(dx * dx + dy * dy);

            auto Nvals = FF(ksi, eta);

            for (int a_index = 0; a_index < edge_nodes.size(); ++a_index) {
                int local_node_index = edge_nodes[a_index];
                double contribution = mult * value * weight * Nvals[local_node_index].real() * detJ_edge;

                F[2 * local_node_index + comp] += contribution;
            }
        }
    }

    return F;
}

template<int NODES>
Eigen::MatrixXd SpectralQuad<NODES>::localC() {
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(nNodes, nNodes);

    for (int j = 0; j < NODES; ++j) {
        double eta = gll_x[j];
        double w_eta = gll_w[j];

        for (int i = 0; i < NODES; ++i) {
            double ksi = gll_x[i];
            double w_ksi = gll_w[i];

            double detJ = std::abs(J(ksi, eta).determinant());
            auto Nvals = FF(ksi, eta);

            for (int a = 0; a < nNodes; ++a)
                for (int b = 0; b < nNodes; ++b)
                    C(a, b) += w_ksi * w_eta * Nvals[a].real() * Nvals[b].real() * detJ;
        }
    }
    return C;
}

template<int NODES>
std::vector<double> SpectralQuad<NODES>::localR(std::vector<double> value) {

    std::vector<double> R(nNodes, 0.0);

    for (int j = 0; j < NODES; ++j) {
        double eta = gll_x[j];
        double w_eta = gll_w[j];

        for (int i = 0; i < NODES; ++i) {
            double ksi = gll_x[i];
            double w_ksi = gll_w[i];

            double detJ = std::abs(J(ksi, eta).determinant());
            auto Nvals = FF(ksi, eta);

            for (int a = 0; a < nNodes; ++a)
                R[a] += value[a] * w_ksi * w_eta * Nvals[a].real() * detJ;
        }
    }
    return R;
}

template<int NODES>
Eigen::MatrixXcd SpectralQuad<NODES>::localM() {
    if (density == 0.0) 
        throw std::runtime_error("Error: density is zero in element " + std::to_string(id));

    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(2 * nNodes, 2 * nNodes);

    for (int j = 0; j < NODES; ++j) {
        double eta = gll_x[j];
        double w_eta = gll_w[j];

        for (int i = 0; i < NODES; ++i) {
            double ksi = gll_x[i];
            double w_ksi = gll_w[i];

            double detJ = std::abs(J(ksi, eta).determinant());
            auto Nvals = FF(ksi, eta);

            for (int a = 0; a < nNodes; ++a)
                for (int b = 0; b < nNodes; ++b) {
                    std::complex<double> Mij = w_ksi * w_eta * Nvals[a] * Nvals[b] * detJ;
                    M(2 * a, 2 * b) += Mij;
                    M(2 * a + 1, 2 * b + 1) += Mij;
                }
        }
    }
    return density * M;
}

template<int NODES>
std::vector<int> SpectralQuad<NODES>::edge_to_node(int edge) {
    std::vector<int> res;

    if (edge < 0 || edge > 3)
        throw std::runtime_error("edge_to_node: wrong edge index");

    if (edge == 0) { // bottom edge: nodes 0, 1, 2, ..., NODES-1
        for (int i = 0; i < NODES; ++i) {
            res.push_back(i);
        }
    }
    else if (edge == 1) { // right edge: nodes NODES-1, 2*NODES-1, 3*NODES-1, ..., NODES*NODES-1
        for (int j = 0; j < NODES; ++j) {
            int node_idx = j * NODES + (NODES - 1);
            res.push_back(node_idx);
        }
    }
    else if (edge == 2) { // top edge: nodes NODES*(NODES-1), NODES*(NODES-1)+1, ..., NODES*NODES-1
        for (int i = 0; i < NODES; ++i) {
            int node_idx = (NODES - 1) * NODES + i;
            res.push_back(node_idx);
        }
    }
    else if (edge == 3) { // left edge: nodes 0, NODES, 2*NODES, ..., NODES*(NODES-1)
        for (int j = 0; j < NODES; ++j) {
            int node_idx = j * NODES;
            res.push_back(node_idx);
        }
    }

    return res;
}
template<int NODES>
double SpectralQuad<NODES>::gaussPoint(LocVar var, int i)
{
    if (i < 0 || i >= NODES * NODES)
        throw std::out_of_range("SpectralQuad<NODES>::gaussPoint: index out of range");

    int idx_eta = i / NODES;
    int idx_ksi = i % NODES;

    if (var == KSI)
        return gll_x[idx_ksi];
    else if (var == ETA)
        return gll_x[idx_eta];
    else
        return 0.0;
}

template<int NODES>
double SpectralQuad<NODES>::weight(LocVar var, int i)
{
    if (i < 0 || i >= NODES * NODES)
        throw std::out_of_range("SpectralQuad<NODES>::weight: index out of range");

    int idx_eta = i / NODES;
    int idx_ksi = i % NODES;

    if (var == KSI)
        return gll_w[idx_ksi];
    else if (var == ETA)
        return gll_w[idx_eta];
    else
        return 0.0;
}

template<int NODES>
double SpectralQuad<NODES>::Volume() {
    std::complex<double> S = 0.0;

    for (int j = 0; j < NODES; ++j) {
        double eta = gll_x[j];
        double w_eta = gll_w[j];

        for (int i = 0; i < NODES; ++i) {
            double ksi = gll_x[i];
            double w_ksi = gll_w[i];

            auto Nvals = FF(ksi, eta);
            double detJ = std::abs(J(ksi, eta).determinant());

            for (int a = 0; a < nNodes; ++a)
                S += w_ksi * w_eta * Nvals[a] * detJ;
        }
    }
    return S.real();
}

template<int NODES>
bool SpectralQuad<NODES>::pointInElem(std::vector<double> point) {
    return false;
}

template<int NODES>
std::vector<double> SpectralQuad<NODES>::coordFF(double x0, double y0, double) {
    return std::vector<double>();
}

template<int NODES>
void SpectralQuad<NODES>::set_pressure(int edge, double value) {
    std::vector<int> nodes_on_edge = edge_to_node(edge);

    if (nodes_on_edge.size() < 2)
        throw std::runtime_error("Not enough nodes on edge");

    int n0 = nodes_on_edge.front();
    int n1 = nodes_on_edge.back();

    double dx = x[n1] - x[n0];
    double dy = y[n1] - y[n0];

    double nx = dy;
    double ny = -dx;

    double length = std::sqrt(nx * nx + ny * ny);
    if (length > 0) {
        nx /= length;
        ny /= length;
    }

    std::pair<int, int> key_x(edge, 0);
    std::pair<int, int> key_y(edge, 1);
    load[key_x] += -value * nx;
    load[key_y] += -value * ny;
}

template<int NODES>
double SpectralQuad<NODES>::len_edge(int edge) {
    auto nodes_on_edge = edge_to_node(edge);

    int n0 = nodes_on_edge.front();
    int n1 = nodes_on_edge.back();

    return std::sqrt(std::pow(x[n0] - x[n1], 2) + std::pow(y[n0] - y[n1], 2));
}