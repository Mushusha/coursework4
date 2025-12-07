#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cassert>

#include "Eigen/Core"
#include "HexSEM.h"
#include "MathMV.h"


template <int NODES>
class SpectralInfHex : public SpectralHex<NODES> {
    static_assert(NODES >= 2, "NODES must be >= 2");

public:
    SpectralInfHex() : SpectralHex<NODES>(), is_dynamic(false), omega(0.0) {
        init_gll();
    }
    SpectralInfHex(int id, ElemType type, std::vector<int> nodes)
        : SpectralHex<NODES>(id, type, nodes), is_dynamic(false), omega(0.0) {
        init_gll();
    }
    SpectralInfHex(const SpectralInfHex& other) 
        : SpectralHex<NODES>(other), is_dynamic(other.is_dynamic), omega(other.omega) {
        init_gll();
    }
    SpectralInfHex& operator=(const SpectralInfHex& other) {
        if (this != &other) {
            SpectralHex<NODES>::operator=(other);
            is_dynamic = other.is_dynamic;
            omega = other.omega;
        }
        return *this;
    }
    SpectralInfHex(SpectralInfHex&& other) noexcept 
        : SpectralHex<NODES>(std::move(other)),
          is_dynamic(other.is_dynamic), omega(other.omega) {
        init_gll();
    }
    SpectralInfHex& operator=(SpectralInfHex&& other) noexcept {
        if (this != &other) {
            SpectralHex<NODES>::operator=(std::move(other));
            is_dynamic = other.is_dynamic;
            omega = other.omega;
        }
        return *this;
    }
    virtual ~SpectralInfHex() = default;

    bool is_dynamic = false;
    double omega = 0.0;
    
    double pole_x = 0.0, pole_y = 0.0, pole_z = 0.0;

    std::vector<std::complex<double>> FF(double ksi, double eta, double zeta) override;
    Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta) override;

    double gaussPoint(LocVar var, int i) override;
    double weight(LocVar var, int i) override;

private:
    static constexpr int nNodes = NODES * NODES * NODES;

    std::vector<double> gr_x;
    std::vector<double> gr_w;

    void init_gll();

    double lagrange1D_ksi(int k, double ksi) const;
    double dlagrange1D_ksi(int k, double ksi) const;

    double decay_function(int k, double ksi) const;
    double ddecay_function(int k, double ksi) const;

    std::complex<double> dynamic_multiplier(double ksi, double eta, double zeta) const;
};


template<int NODES>
void SpectralInfHex<NODES>::init_gll() {
    SpectralHex<NODES>::init_gll();

    compute_gauss_radau_nodes_weights(NODES, gr_x, gr_w);
}

template<int NODES>
double SpectralInfHex<NODES>::lagrange1D_ksi(int k, double ksi) const {
    double val = 1.0;
    for (int m = 0; m < NODES; ++m) {
        if (m == k) continue;
        val *= (ksi - gr_x[m]) / (gr_x[k] - gr_x[m]);
    }
    return val;
}

template<int NODES>
double SpectralInfHex<NODES>::dlagrange1D_ksi(int k, double ksi) const {
    double sum = 0.0;
    for (int m = 0; m < NODES; ++m) {
        if (m == k) continue;
        double prod = 1.0 / (gr_x[k] - gr_x[m]);
        for (int n = 0; n < NODES; ++n) {
            if (n == k || n == m) continue;
            prod *= (ksi - gr_x[n]) / (gr_x[k] - gr_x[n]);
        }
        sum += prod;
    }
    return sum;
}

template<int NODES>
double SpectralInfHex<NODES>::decay_function(int k, double ksi) const {
    double L_k = lagrange1D_ksi(k, ksi);
    double decay = 2.0 / (1.0 - ksi);
    return L_k * decay;
}

template<int NODES>
double SpectralInfHex<NODES>::ddecay_function(int k, double ksi) const {
    double L_k = lagrange1D_ksi(k, ksi);
    double dL_k = dlagrange1D_ksi(k, ksi);

    double one_minus_ksi = 1.0 - ksi;
    double decay = 2.0 / one_minus_ksi;
    double d_decay = 2.0 / (one_minus_ksi * one_minus_ksi);

    return dL_k * decay + L_k * d_decay;
}

template<int NODES>
std::complex<double> SpectralInfHex<NODES>::dynamic_multiplier(
    double ksi, double eta, double zeta) const {
    if (!is_dynamic || omega <= 0.0) {
        return std::complex<double>(1.0, 0.0);
    }

    double lambda = this->Young * this->Poisson / ((1.0 + this->Poisson) * (1.0 - 2.0 * this->Poisson));
    double mu = this->Young / (2.0 * (1.0 + this->Poisson));
    double c = std::sqrt((lambda + 2.0 * mu) / this->density);
    
    double k_wave = omega / c;
    
    double cx = 0.0, cy = 0.0, cz = 0.0;
    int face_nodes = NODES * NODES;
    for (int i = 0; i < face_nodes; ++i) {
        cx += this->x[i];
        cy += this->y[i];
        cz += this->z[i];
    }
    cx /= face_nodes; cy /= face_nodes; cz /= face_nodes;
    double A_len = std::sqrt((cx - pole_x)*(cx - pole_x) + 
                             (cy - pole_y)*(cy - pole_y) + 
                             (cz - pole_z)*(cz - pole_z));
    if (A_len < 1e-10) A_len = 1.0;
    
    double one_minus_ksi = 1.0 - ksi;
    std::complex<double> i(0.0, 1.0);
    double decay_factor = std::sqrt(2.0 / one_minus_ksi);
    std::complex<double> phase1 = std::exp(i * k_wave * A_len / 2.0);
    std::complex<double> phase2 = std::exp(i * k_wave * A_len / one_minus_ksi);
    
    return decay_factor * phase1 * phase2;
}

template<int NODES>
std::vector<std::complex<double>> SpectralInfHex<NODES>::FF(
    double ksi, double eta, double zeta) {
    std::vector<std::complex<double>> N(nNodes);
    int idx = 0;

    std::complex<double> dyn_mult = dynamic_multiplier(ksi, eta, zeta);

    for (int k = 0; k < NODES; ++k) {
        double L_k = this->lagrange1D(k, zeta);
        for (int j = 0; j < NODES; ++j) {
            double L_j = this->lagrange1D(j, eta);
            for (int i = 0; i < NODES; ++i) {
                double M_i = decay_function(i, ksi);
                N[idx] = dyn_mult * M_i * L_j * L_k;
                idx++;
            }
        }
    }
    return N;
}

template<int NODES>
Eigen::MatrixXcd SpectralInfHex<NODES>::gradFF(
    double ksi, double eta, double zeta) {
    Eigen::MatrixXcd grad(3, nNodes);
    grad.setZero();
    int idx = 0;

    std::complex<double> dyn_mult = dynamic_multiplier(ksi, eta, zeta);

    for (int k = 0; k < NODES; ++k) {
        double L_k = this->lagrange1D(k, zeta);
        double dL_k = this->dlagrange1D(k, zeta);
        for (int j = 0; j < NODES; ++j) {
            double L_j = this->lagrange1D(j, eta);
            double dL_j = this->dlagrange1D(j, eta);
            for (int i = 0; i < NODES; ++i, ++idx) {
                double M_i = decay_function(i, ksi);
                double dM_i = ddecay_function(i, ksi);

                grad(0, idx) = dyn_mult * dM_i * L_j * L_k;
                grad(1, idx) = dyn_mult * M_i * dL_j * L_k;
                grad(2, idx) = dyn_mult * M_i * L_j * dL_k;
            }
        }
    }
    return grad;
}

template<int NODES>
double SpectralInfHex<NODES>::gaussPoint(LocVar var, int i) {
    if (i < 0 || i >= nNodes)
        throw std::out_of_range("SpectralInfHex::gaussPoint: index out of range");

    int idx_ksi = i % NODES;
    int idx_eta = (i / NODES) % NODES;
    int idx_zeta = i / (NODES * NODES);

    if (var == KSI)
        return gr_x[idx_ksi];
    else if (var == ETA)
        return this->gll_x[idx_eta];
    else if (var == ZETA)
        return this->gll_x[idx_zeta];
    else
        return 0.0;
}

template<int NODES>
double SpectralInfHex<NODES>::weight(LocVar var, int i) {
    if (i < 0 || i >= nNodes)
        throw std::out_of_range("SpectralInfHex::weight: index out of range");

    int idx_ksi = i % NODES;
    int idx_eta = (i / NODES) % NODES;
    int idx_zeta = i / (NODES * NODES);

    if (var == KSI)
        return gr_w[idx_ksi];
    else if (var == ETA)
        return this->gll_w[idx_eta];
    else if (var == ZETA)
        return this->gll_w[idx_zeta];
    else
        return 0.0;
}