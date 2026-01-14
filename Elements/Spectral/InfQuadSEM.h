#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

#include "Eigen/Core"
#include "QuadSEM.h"
#include "MathMV.h"
#include "Element.h"

class Data;

namespace {
    inline std::string format_k_matrix_stats_inf(const Eigen::MatrixXcd& K, int elem_id, const char* elem_name) {
        const int rows = static_cast<int>(K.rows());
        const int cols = static_cast<int>(K.cols());
        const int total = rows * cols;

        int non_zero = 0;
        double max_abs = 0.0;
        double min_abs = std::numeric_limits<double>::infinity();
        double sum_abs = 0.0;

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                const double a = std::abs(K(r, c));
                if (a != 0.0) ++non_zero;
                if (a > max_abs) max_abs = a;
                if (a < min_abs) min_abs = a;
                sum_abs += a;
            }
        }
        if (!std::isfinite(min_abs)) min_abs = 0.0;

        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(6);
        oss << elem_name << " element ID=" << elem_id
            << " K matrix: size=" << rows << "x" << cols
            << ", non_zero=" << non_zero
            << ", max_abs=" << max_abs
            << ", min_abs=" << min_abs
            << ", avg_abs=" << (total > 0 ? (sum_abs / static_cast<double>(total)) : 0.0);
        return oss.str();
    }
}


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
    Eigen::MatrixXcd J(double ksi, double eta, double zeta = 0) override;
    Eigen::MatrixXcd B(double ksi = 0, double eta = 0, double zeta = 0) override;
    Eigen::MatrixXcd localK() override;
    Eigen::MatrixXd localDamping() override;

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
    std::complex<double> ddynamic_multiplier_dksi(double ksi, double eta) const;
    
    double compute_A_len() const;
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
    
    double one_minus_ksi = 1.0 - ksi;
    if (one_minus_ksi < 1e-10) one_minus_ksi = 1e-10;
    
    double decay;
    decay = 1.0 / one_minus_ksi;
    
    /*
    if (Data::is_dynamic && omega > 0.0) {
        double cx = 0.0, cy = 0.0;
        for (int k = 0; k < NODES; ++k) {
            int idx = 0 + k * NODES;
            cx += this->x[idx];
            cy += this->y[idx];
        }
        cx /= NODES;
        cy /= NODES;
        double A_len = std::sqrt((cx - pole_x) * (cx - pole_x) + (cy - pole_y) * (cy - pole_y));
        if (A_len < 1e-10) A_len = 1.0;
        
        decay = std::sqrt(2.0 * A_len) / std::sqrt(one_minus_ksi);
    } else {
        decay = 1.0 / one_minus_ksi;
    }
    */
    
    return L_j * decay;
}

template<int NODES>
double SpectralInfQuad<NODES>::ddecay_function(int j, double ksi) const {
    double L_j = lagrange1D_ksi(j, ksi);
    double dL_j = dlagrange1D_ksi(j, ksi);
    
    double one_minus_ksi = 1.0 - ksi;
    if (one_minus_ksi < 1e-10) one_minus_ksi = 1e-10;
    
    double decay, d_decay;
    decay = 1.0 / one_minus_ksi;
    d_decay = 1.0 / (one_minus_ksi * one_minus_ksi);
    
    /* 
    if (Data::is_dynamic && omega > 0.0) {
        double cx = 0.0, cy = 0.0;
        for (int k = 0; k < NODES; ++k) {
            int idx = 0 + k * NODES;
            cx += this->x[idx];
            cy += this->y[idx];
        }
        cx /= NODES;
        cy /= NODES;
        double A_len = std::sqrt((cx - pole_x) * (cx - pole_x) + (cy - pole_y) * (cy - pole_y));
        if (A_len < 1e-10) A_len = 1.0;
        
        decay = std::sqrt(2.0 * A_len) / std::sqrt(one_minus_ksi);
        d_decay = 0.5 * std::sqrt(2.0 * A_len) / std::pow(one_minus_ksi, 1.5);
    } else {
        decay = 1.0 / one_minus_ksi;
        d_decay = 1.0 / (one_minus_ksi * one_minus_ksi);
    }
    */
    
    return dL_j * decay + L_j * d_decay;
}

template<int NODES>
std::complex<double> SpectralInfQuad<NODES>::dynamic_multiplier(double ksi, double eta) const {
    return std::complex<double>(1.0, 0.0);
    
    /* 
    if (!is_dynamic || omega <= 0.0) {
        return std::complex<double>(1.0, 0.0);
    }
    
    double lambda = this->Young * this->Poisson / ((1.0 + this->Poisson) * (1.0 - 2.0 * this->Poisson));
    double mu = this->Young / (2.0 * (1.0 + this->Poisson));
    double c = std::sqrt((lambda + 2.0 * mu) / this->density);
    double k_wave = omega / c;
    
    double cx = 0.0, cy = 0.0;
    for (int k = 0; k < NODES; ++k) {
        int idx = 0 + k * NODES;
        cx += this->x[idx];
        cy += this->y[idx];
    }
    cx /= NODES;
    cy /= NODES;
    double A_len = std::sqrt((cx - pole_x) * (cx - pole_x) + (cy - pole_y) * (cy - pole_y));
    if (A_len < 1e-10) A_len = 1.0;
    
    double one_minus_ksi = 1.0 - ksi;
    if (one_minus_ksi < 1e-10) one_minus_ksi = 1e-10;
    double r = A_len / one_minus_ksi;
    std::complex<double> phase = std::exp(std::complex<double>(0.0, 1.0) * k_wave * A_len * 0.5) *
                                  std::exp(std::complex<double>(0.0, 1.0) * k_wave * r);
    
    return phase;
    */
}

template<int NODES>
std::complex<double> SpectralInfQuad<NODES>::ddynamic_multiplier_dksi(double ksi, double eta) const {
    return std::complex<double>(0.0, 0.0);
    
    /* 
    if (!is_dynamic || omega <= 0.0) {
        return std::complex<double>(0.0, 0.0);
    }
    
    double lambda = this->Young * this->Poisson / ((1.0 + this->Poisson) * (1.0 - 2.0 * this->Poisson));
    double mu = this->Young / (2.0 * (1.0 + this->Poisson));
    double c = std::sqrt((lambda + 2.0 * mu) / this->density);
    double k_wave = omega / c;
    
    double cx = 0.0, cy = 0.0;
    for (int k = 0; k < NODES; ++k) {
        int idx = 0 + k * NODES;
        cx += this->x[idx];
        cy += this->y[idx];
    }
    cx /= NODES;
    cy /= NODES;
    double A_len = std::sqrt((cx - pole_x) * (cx - pole_x) + (cy - pole_y) * (cy - pole_y));
    if (A_len < 1e-10) A_len = 1.0;
    
    double one_minus_ksi = 1.0 - ksi;
    if (one_minus_ksi < 1e-10) one_minus_ksi = 1e-10;
    
    std::complex<double> phase = std::exp(std::complex<double>(0.0, 1.0) * k_wave * A_len * 0.5) *
                                  std::exp(std::complex<double>(0.0, 1.0) * k_wave * A_len / one_minus_ksi);
    std::complex<double> dphase = phase * std::complex<double>(0.0, 1.0) * k_wave * A_len / (one_minus_ksi * one_minus_ksi);
    
    return dphase;
    */
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
    std::complex<double> ddyn_mult_dksi = ddynamic_multiplier_dksi(ksi, eta);

    for (int j = 0; j < NODES; ++j) {
        double L_j = this->lagrange1D(j, eta);
        double dL_j = this->dlagrange1D(j, eta);
        for (int i = 0; i < NODES; ++i, ++idx) {
            double M_i = decay_function(i, ksi);
            double dM_i = ddecay_function(i, ksi);

            grad(0, idx) = ddyn_mult_dksi * M_i * L_j + dyn_mult * dM_i * L_j;
            grad(1, idx) = dyn_mult * M_i * dL_j;
        }
    }
    return grad;
}

template<int NODES>
Eigen::MatrixXcd SpectralInfQuad<NODES>::J(double ksi, double eta, double) {
    Eigen::MatrixXcd jac = Eigen::MatrixXcd::Zero(2, 2);
    
    int idx = 0;
    for (int j = 0; j < NODES; ++j) {
        double L_j = this->lagrange1D(j, eta);
        double dL_j = this->dlagrange1D(j, eta);
        for (int i = 0; i < NODES; ++i, ++idx) {
            double M_i = decay_function(i, ksi);
            double dM_i = ddecay_function(i, ksi);
            
            double dN_dksi = dM_i * L_j;
            double dN_deta = M_i * dL_j;
            
            jac(0, 0) += std::complex<double>(dN_dksi, 0.0) * this->x[idx];
            jac(0, 1) += std::complex<double>(dN_dksi, 0.0) * this->y[idx];
            jac(1, 0) += std::complex<double>(dN_deta, 0.0) * this->x[idx];
            jac(1, 1) += std::complex<double>(dN_deta, 0.0) * this->y[idx];
        }
    }
    return jac;
}

template<int NODES>
Eigen::MatrixXcd SpectralInfQuad<NODES>::B(double ksi, double eta, double) {
    Eigen::MatrixXcd Bm = Eigen::MatrixXcd::Zero(3, 2 * nNodes);
    Eigen::Matrix2cd invJ = this->J(ksi, eta).inverse();
    Eigen::MatrixXcd dN = invJ * gradFF(ksi, eta);

    for (int a = 0; a < nNodes; ++a) {
        std::complex<double> dNx = dN(X, a);
        std::complex<double> dNy = dN(Y, a);

        Bm(0, 2 * a) = dNx;
        Bm(1, 2 * a + 1) = dNy;
        Bm(2, 2 * a) = dNy;
        Bm(2, 2 * a + 1) = dNx;
    }
    return Bm;
}

template<int NODES>
Eigen::MatrixXcd SpectralInfQuad<NODES>::localK() {
    Eigen::MatrixXcd K = Eigen::MatrixXcd::Zero(2 * nNodes, 2 * nNodes);
    
    for (int j = 0; j < NODES; ++j) {
        for (int i = 0; i < NODES; ++i) {
            int node_idx = i + j * NODES;
            double ksi = gaussPoint(KSI, node_idx);
            double eta = gaussPoint(ETA, node_idx);
            double w_ksi = weight(KSI, node_idx);
            double w_eta = weight(ETA, node_idx);

            Eigen::MatrixXcd Bm = B(ksi, eta);
            Eigen::MatrixXcd Jm = this->J(ksi, eta);

            double detJ = std::abs(Jm.determinant());

            Eigen::MatrixXcd BtDB = Bm.adjoint() * this->D.template cast<std::complex<double>>() * Bm;
            Eigen::MatrixXcd contribution = (w_ksi * w_eta) * BtDB * detJ;
            
            K += contribution;
        }
    }
    
    double A_len = compute_A_len();
    if (A_len < 1e-10) A_len = 1.0;
    
    double reference_size = Data::reference_element_size; 
    if (reference_size < 1e-10) reference_size = 1.0;
    
    const double base_scale = 1e-3;
    double scale_factor = base_scale * (reference_size / A_len);
    
    if (scale_factor > 1.0) scale_factor = 1.0;
    if (scale_factor < 1e-5) scale_factor = 1e-5;
    
    K = K * scale_factor;
    
    return K;
}

template<int NODES>
Eigen::MatrixXd SpectralInfQuad<NODES>::localDamping() {
    Eigen::MatrixXd Dmp = Eigen::MatrixXd::Zero(nNodes, nNodes);

    if (!Data::is_dynamic || omega <= 0.0) {
        return Dmp;
    }

    const double nu = this->Poisson;
    const double E = this->Young;
    const double rho = this->density;
    if (rho <= 0.0) return Dmp;

    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));
    const double c_p = std::sqrt((lambda + 2.0 * mu) / rho);
    const double c_s = std::sqrt(mu / rho);

    const double Z_p = rho * c_p;
    const double Z_s = rho * c_s;
    const double Z_eff = 0.75 * Z_p + 0.25 * Z_s;
    
    double frequency_factor = 1.0;
    if (omega > 0.0) {
        if (omega < 30.0) {
            frequency_factor = 1.0 + 2.0 * (30.0 - omega) / 30.0;
        } else if (omega < 100.0) {
            frequency_factor = 1.0 + (100.0 - omega) / 100.0;
        } else {
            frequency_factor = 0.8;
        }
    }
    
    const double base_damping_factor = 4.5;
    const double damping_factor = base_damping_factor * frequency_factor;
    const double c_eff = Z_eff * damping_factor;

    const double ksi = -1.0;

    double total_damping = 0.0;
    int integration_points = 0;

    for (int j = 0; j < NODES; ++j) {
        const double eta = this->gll_x[j];
        const double w_eta = this->gll_w[j];

        Eigen::MatrixXcd grad = gradFF(ksi, eta);
        double dx_eta = 0.0, dy_eta = 0.0;
        for (int a = 0; a < nNodes; ++a) {
            dx_eta += grad(ETA, a).real() * this->x[a];
            dy_eta += grad(ETA, a).real() * this->y[a];
        }
        const double detJ_edge = std::sqrt(dx_eta * dx_eta + dy_eta * dy_eta);

        if (detJ_edge < 1e-10) continue;

        double nx = -dy_eta / detJ_edge;
        double ny = dx_eta / detJ_edge;
        
        double cx = 0.0, cy = 0.0;
        int boundary_node_count = 0;
        for (int a = 0; a < nNodes; ++a) {
            int i_a = a % NODES;
            if (i_a == 0) {
                cx += this->x[a];
                cy += this->y[a];
                boundary_node_count++;
            }
        }
        if (boundary_node_count > 0) {
            cx /= boundary_node_count;
            cy /= boundary_node_count;
            
            double dx_pole = cx - pole_x;
            double dy_pole = cy - pole_y;
            
            double dot = nx * dx_pole + ny * dy_pole;
            if (dot < 0.0) {
                nx = -nx;
                ny = -ny;
            }
        }

        auto Nvals = FF(ksi, eta);

        for (int a = 0; a < nNodes; ++a) {
            int i_a = a % NODES;
            bool is_boundary_node_a = (i_a == 0);
            
            if (!is_boundary_node_a) continue;
            
            const double Na = Nvals[a].real();
            
            for (int b = 0; b < nNodes; ++b) {
                int i_b = b % NODES;
                bool is_boundary_node_b = (i_b == 0);
                
                if (!is_boundary_node_b) continue;
                
                const double Nb = Nvals[b].real();
                
                const double contribution = w_eta * Na * Nb * detJ_edge * c_eff;
                Dmp(a, b) += contribution;
                total_damping += std::abs(contribution);
            }
        }
        integration_points++;
    }

    return Dmp;
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

template<int NODES>
double SpectralInfQuad<NODES>::compute_A_len() const {
    double cx = 0.0, cy = 0.0;
    for (int k = 0; k < NODES; ++k) {
        int idx = 0 + k * NODES;
        cx += this->x[idx];
        cy += this->y[idx];
    }
    cx /= NODES;
    cy /= NODES;
    double A_len = std::sqrt((cx - pole_x) * (cx - pole_x) + (cy - pole_y) * (cy - pole_y));
    return A_len;
}

