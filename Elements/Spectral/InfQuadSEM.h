#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cassert>

#include "Eigen/Core"
#include "QuadSEM.h"
#include "MathMV.h"


template <int NODES>
class SpectralInfQuad : public SpectralQuad<NODES> {
    static_assert(NODES >= 2, "NODES must be >= 2");

public:
    SpectralInfQuad() : SpectralQuad<NODES>(), is_dynamic(false), omega(0.0) {
        init_gll();
    }
    SpectralInfQuad(int id, ElemType type, std::vector<int> nodes)
        : SpectralQuad<NODES>(id, type, nodes), is_dynamic(false), omega(0.0) {
        init_gll();
    }
    SpectralInfQuad(const SpectralInfQuad& other) 
        : SpectralQuad<NODES>(other), is_dynamic(other.is_dynamic), omega(other.omega) {
        init_gll();
    }
    SpectralInfQuad& operator=(const SpectralInfQuad& other) {
        if (this != &other) {
            SpectralQuad<NODES>::operator=(other);
            is_dynamic = other.is_dynamic;
            omega = other.omega;
        }
        return *this;
    }
    SpectralInfQuad(SpectralInfQuad&& other) noexcept 
        : SpectralQuad<NODES>(std::move(other)),
          is_dynamic(other.is_dynamic), omega(other.omega) {
        init_gll();
    }
    SpectralInfQuad& operator=(SpectralInfQuad&& other) noexcept {
        if (this != &other) {
            SpectralQuad<NODES>::operator=(std::move(other));
            is_dynamic = other.is_dynamic;
            omega = other.omega;
        }
        return *this;
    }
    virtual ~SpectralInfQuad() = default;

    bool is_dynamic = false;
    double omega = 0.0;
    
    double pole_x = 0.0, pole_y = 0.0, pole_z = 0.0;

    std::vector<std::complex<double>> FF(double ksi, double eta, double zeta = 0) override;
    Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta = 0) override;

    double gaussPoint(LocVar var, int i) override;
    double weight(LocVar var, int i) override;

private:
    static constexpr int nNodes = NODES * NODES;
    
    std::vector<double> gr_x;
    std::vector<double> gr_w;

    void init_gll();

    double decay_function(int j, double ksi) const;
    double ddecay_function(int j, double ksi) const;

    double lagrange1D_ksi(int j, double ksi) const;
    double dlagrange1D_ksi(int j, double ksi) const;

    std::complex<double> dynamic_multiplier(double ksi, double eta) const;
};


template<int NODES>
void SpectralInfQuad<NODES>::init_gll() {
    SpectralQuad<NODES>::init_gll();

    compute_gauss_radau_nodes_weights(NODES, gr_x, gr_w);
}

template<int NODES>
double SpectralInfQuad<NODES>::lagrange1D_ksi(int j, double ksi) const {
    double val = 1.0;
    for (int m = 0; m < NODES; ++m) {
        if (m == j) continue;
        val *= (ksi - gr_x[m]) / (gr_x[j] - gr_x[m]);
    }
    return val;
}

template<int NODES>
double SpectralInfQuad<NODES>::dlagrange1D_ksi(int j, double ksi) const {
    double sum = 0.0;
    for (int m = 0; m < NODES; ++m) {
        if (m == j) continue;
        double prod = 1.0 / (gr_x[j] - gr_x[m]);
        for (int k = 0; k < NODES; ++k) {
            if (k == j || k == m) continue;
            prod *= (ksi - gr_x[k]) / (gr_x[j] - gr_x[k]);
        }
        sum += prod;
    }
    return sum;
}

template<int NODES>
double SpectralInfQuad<NODES>::decay_function(int j, double ksi) const {
    double L_j = lagrange1D_ksi(j, ksi);
    double decay = 2.0 / (1.0 - ksi);
    return L_j * decay;
}

template<int NODES>
double SpectralInfQuad<NODES>::ddecay_function(int j, double ksi) const {
    double L_j = lagrange1D_ksi(j, ksi);
    double dL_j = dlagrange1D_ksi(j, ksi);
    
    double one_minus_ksi = 1.0 - ksi;
    double decay = 2.0 / one_minus_ksi;
    double d_decay = 2.0 / (one_minus_ksi * one_minus_ksi);
    
    return dL_j * decay + L_j * d_decay;
}

template<int NODES>
std::complex<double> SpectralInfQuad<NODES>::dynamic_multiplier(double ksi, double eta) const {
    if (!is_dynamic || omega <= 0.0) {
        return std::complex<double>(1.0, 0.0);
    }
    
    double lambda = this->Young * this->Poisson / ((1.0 + this->Poisson) * (1.0 - 2.0 * this->Poisson));
    double mu = this->Young / (2.0 * (1.0 + this->Poisson));
    double c = std::sqrt((lambda + 2.0 * mu) / this->density);
    
    double k_wave = omega / c;
    
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (int j = 0; j < NODES; ++j) {
        int idx = j * NODES;
        cx += this->x[idx];
        cy += this->y[idx];
        cz += this->z[idx];
    }
    cx /= NODES; cy /= NODES; cz /= NODES;
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
std::vector<std::complex<double>> SpectralInfQuad<NODES>::FF(double ksi, double eta, double) {
    std::vector<std::complex<double>> N(nNodes);
    int idx = 0;
    
    std::complex<double> dyn_mult = dynamic_multiplier(ksi, eta);

    for (int j = 0; j < NODES; ++j) {
        double L_j = this->lagrange1D(j, eta);
        for (int i = 0; i < NODES; ++i) {
            double M_i = decay_function(i, ksi);
            N[idx] = dyn_mult * M_i * L_j;
            idx++;
        }
    }
    return N;
}

template<int NODES>
Eigen::MatrixXcd SpectralInfQuad<NODES>::gradFF(double ksi, double eta, double) {
    Eigen::MatrixXcd grad(2, nNodes);
    grad.setZero();
    int idx = 0;

    std::complex<double> dyn_mult = dynamic_multiplier(ksi, eta);

    for (int j = 0; j < NODES; ++j) {
        double L_j = this->lagrange1D(j, eta);
        double dL_j = this->dlagrange1D(j, eta);
        for (int i = 0; i < NODES; ++i, ++idx) {
            double M_i = decay_function(i, ksi);
            double dM_i = ddecay_function(i, ksi);

            grad(0, idx) = dyn_mult * dM_i * L_j;
            grad(1, idx) = dyn_mult * M_i * dL_j;
        }
    }
    return grad;
}

template<int NODES>
double SpectralInfQuad<NODES>::gaussPoint(LocVar var, int i) {
    if (i < 0 || i >= NODES * NODES)
        throw std::out_of_range("SpectralInfQuad::gaussPoint: index out of range");

    int idx_eta = i / NODES;
    int idx_ksi = i % NODES;

    if (var == KSI)
        return gr_x[idx_ksi];
    else if (var == ETA)
        return this->gll_x[idx_eta];
    else
        return 0.0;
}

template<int NODES>
double SpectralInfQuad<NODES>::weight(LocVar var, int i) {
    if (i < 0 || i >= NODES * NODES)
        throw std::out_of_range("SpectralInfQuad::weight: index out of range");

    int idx_eta = i / NODES;
    int idx_ksi = i % NODES;

    if (var == KSI)
        return gr_w[idx_ksi];
    else if (var == ETA)
        return this->gll_w[idx_eta];
    else
        return 0.0;
}
