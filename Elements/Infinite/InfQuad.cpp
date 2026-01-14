#include "InfQuad.h"

double InfQuad::compute_A() const {
	double cx = (x[0] + x[3]) / 2.0;
	double cy = (y[0] + y[3]) / 2.0;
	double cz = (z[0] + z[3]) / 2.0;
	
	double dx = cx - pole_x;
	double dy = cy - pole_y;
	double dz = cz - pole_z;
	
	return std::sqrt(dx*dx + dy*dy + dz*dz);
}

std::vector<std::complex<double>> InfQuad::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> N;
	N.resize(4);

	double one_minus_ksi = 1.0 - ksi;
	if (one_minus_ksi < 1e-10) one_minus_ksi = 1e-10;
	
	double M_Q, M_C;
	std::complex<double> phase_mult(1.0, 0.0);
	
	M_Q = 1.0 / one_minus_ksi; 
	M_C = -ksi / one_minus_ksi;
	
	/* 
	if (is_dynamic && omega > 0.0) {
		double A_len = compute_A();
		if (A_len < 1e-10) A_len = 1.0;
		
		double decay_factor = std::sqrt(2.0 * A_len) / std::sqrt(one_minus_ksi);
		
		M_Q = decay_factor / one_minus_ksi;
		M_C = -ksi * decay_factor / one_minus_ksi;
		
		double lambda = get_E() * get_nu() / ((1.0 + get_nu()) * (1.0 - 2.0 * get_nu()));
		double mu = get_E() / (2.0 * (1.0 + get_nu()));
		double c = std::sqrt((lambda + 2.0 * mu) / get_rho());
		double k_wave = omega / c;
		double r = A_len / one_minus_ksi;
		phase_mult = std::exp(std::complex<double>(0.0, 1.0) * k_wave * A_len * 0.5) *
		              std::exp(std::complex<double>(0.0, 1.0) * k_wave * r);
	} else {
		M_Q = 1.0 / one_minus_ksi; 
		M_C = -ksi / one_minus_ksi;
	}
	*/
	
	double L_eta_minus = (1.0 - eta) / 2.0;
	double L_eta_plus = (1.0 + eta) / 2.0;
	
	N[0] = phase_mult * M_Q * L_eta_minus;
	N[1] = phase_mult * M_C * L_eta_minus;
	N[2] = phase_mult * M_C * L_eta_plus;
	N[3] = phase_mult * M_Q * L_eta_plus;

	return N;
}

Eigen::MatrixXcd InfQuad::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd grad = Eigen::MatrixXcd::Zero(2, 4);
	
	double one_minus_ksi = 1.0 - ksi;
	if (one_minus_ksi < 1e-10) one_minus_ksi = 1e-10;
	
	double M_Q, M_C;
	double dM_Q_dksi, dM_C_dksi;
	std::complex<double> phase_mult(1.0, 0.0);
	std::complex<double> dphase_mult_dksi(0.0, 0.0);
	
	double d = one_minus_ksi * one_minus_ksi;
	M_Q = 1.0 / one_minus_ksi;
	M_C = -ksi / one_minus_ksi;
	dM_Q_dksi = 1.0 / d;
	dM_C_dksi = -(1.0 + ksi) / d;
	
	/*
	if (is_dynamic && omega > 0.0) {
		double A_len = compute_A();
		if (A_len < 1e-10) A_len = 1.0;
		
		double decay_factor = std::sqrt(2.0 * A_len) / std::sqrt(one_minus_ksi);
		double ddecay_factor = 0.5 * std::sqrt(2.0 * A_len) / std::pow(one_minus_ksi, 1.5);
		
		double d = one_minus_ksi * one_minus_ksi;

		M_Q = decay_factor / one_minus_ksi;
		M_C = -ksi * decay_factor / one_minus_ksi;
		dM_Q_dksi = ddecay_factor / one_minus_ksi + decay_factor / d;
		dM_C_dksi = -decay_factor / one_minus_ksi - ksi * ddecay_factor / one_minus_ksi - ksi * decay_factor / d;
		
		double lambda = get_E() * get_nu() / ((1.0 + get_nu()) * (1.0 - 2.0 * get_nu()));
		double mu = get_E() / (2.0 * (1.0 + get_nu()));
		double c = std::sqrt((lambda + 2.0 * mu) / get_rho());
		double k_wave = omega / c;
		double r = A_len / one_minus_ksi;
		phase_mult = std::exp(std::complex<double>(0.0, 1.0) * k_wave * A_len * 0.5) *
		              std::exp(std::complex<double>(0.0, 1.0) * k_wave * r);
		dphase_mult_dksi = phase_mult * std::complex<double>(0.0, 1.0) * k_wave * A_len / d;
	} else {
		double d = one_minus_ksi * one_minus_ksi;
		M_Q = 1.0 / one_minus_ksi;
		M_C = -ksi / one_minus_ksi;
		dM_Q_dksi = 1.0 / d;
		dM_C_dksi = -(1.0 + ksi) / d;
	}
	*/
	
	double L_eta_minus = (1.0 - eta) / 2.0;
	double L_eta_plus = (1.0 + eta) / 2.0;
	
	grad(KSI, 0) = (dphase_mult_dksi * M_Q + phase_mult * dM_Q_dksi) * L_eta_minus;
	grad(KSI, 1) = (dphase_mult_dksi * M_C + phase_mult * dM_C_dksi) * L_eta_minus;
	grad(KSI, 2) = (dphase_mult_dksi * M_C + phase_mult * dM_C_dksi) * L_eta_plus;
	grad(KSI, 3) = (dphase_mult_dksi * M_Q + phase_mult * dM_Q_dksi) * L_eta_plus;
	
	grad(ETA, 0) = -phase_mult * M_Q * 0.5;
	grad(ETA, 1) = phase_mult * M_C * 0.5;
	grad(ETA, 2) = -phase_mult * M_C * 0.5;
	grad(ETA, 3) = phase_mult * M_Q * 0.5;
	
	return grad;
}

double InfQuad::gaussPoint(LocVar var, int i) {
	static const double gr_0 = -1.0;
	static const double gr_1 = 1.0 / 3.0;
	
	static const double gp_std = 1.0 / std::sqrt(3.0);
	
	int eta_idx = i / 2;
	int ksi_idx = i % 2;
	
	if (var == KSI) {
		return (ksi_idx == 0) ? gr_1 : gr_0;
	}
	else {
		return (eta_idx == 0) ? -gp_std : gp_std;
	}
}

double InfQuad::weight(LocVar var, int i) {
	
	int ksi_idx = i % 2;
	
	if (var == KSI) {  
		return (ksi_idx == 0) ? 1.5 : 0.5;
	}
	else {
		return 1.0;
	}
}

