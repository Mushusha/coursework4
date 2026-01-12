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
	double M_Q = 1.0 / one_minus_ksi; 
	double M_C = -ksi / one_minus_ksi;
	
	double L_eta_minus = (1.0 - eta) / 2.0;
	double L_eta_plus = (1.0 + eta) / 2.0;
	
	N[0] = M_Q * L_eta_minus;
	N[1] = M_C * L_eta_minus;
	N[2] = M_C * L_eta_plus;
	N[3] = M_Q * L_eta_plus;

	if (is_dyn && omega > 0.0) {
		double lambda = Young * Poisson / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
		double mu = Young / (2.0 * (1.0 + Poisson));
		double c = std::sqrt((lambda + 2.0 * mu) / density);
		
		double k_wave = omega / c;
		
		double A_len = compute_A();
		
		std::complex<double> i(0.0, 1.0);

		double decay_factor = std::sqrt(2.0 / one_minus_ksi);

		std::complex<double> phase1 = std::exp(i * k_wave * A_len / 2.0);
		std::complex<double> phase2 = std::exp(i * k_wave * A_len * ksi / 2.0);
		std::complex<double> mult = decay_factor * phase1 * phase2;

		std::transform(N.begin(), N.end(), N.begin(), 
			[&mult](std::complex<double> a) { return a * mult; });
	}

	return N;
}

Eigen::MatrixXcd InfQuad::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd grad = Eigen::MatrixXcd::Zero(2, 4);
	
	double one_minus_ksi = 1.0 - ksi;
	double d = one_minus_ksi * one_minus_ksi;
	
	double dM_Q_dksi = 1.0 / d;
	double dM_C_dksi = -(1.0 + ksi) / d;
	
	grad(KSI, 0) =  dM_Q_dksi * (1.0 - eta) / 2.0;
	grad(KSI, 1) =  dM_C_dksi * (1.0 - eta) / 2.0;
	grad(KSI, 2) =  dM_C_dksi * (1.0 + eta) / 2.0;
	grad(KSI, 3) =  dM_Q_dksi * (1.0 + eta) / 2.0;
	
	grad(ETA, 0) = -0.5 / one_minus_ksi;
	grad(ETA, 1) =  0.5 * ksi / one_minus_ksi;
	grad(ETA, 2) = -0.5 * ksi / one_minus_ksi;
	grad(ETA, 3) =  0.5 / one_minus_ksi;
	
	if (is_dyn && omega > 0.0) {
		double lambda = Young * Poisson / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
		double mu = Young / (2.0 * (1.0 + Poisson));
		double c = std::sqrt((lambda + 2.0 * mu) / density);
		double k_wave = omega / c;
		double A_len = compute_A();
		
		std::complex<double> i(0.0, 1.0);
		double decay_factor = std::sqrt(2.0 / one_minus_ksi);
		double ddecay_factor = 0.5 * std::sqrt(2.0) * std::pow(one_minus_ksi, -1.5);
		
		std::complex<double> phase1 = std::exp(i * k_wave * A_len / 2.0);
		std::complex<double> phase2 = std::exp(i * k_wave * A_len * ksi / 2.0);
		std::complex<double> dphase2 = phase2 * i * k_wave * A_len / 2.0;
		
		std::complex<double> dyn_mult = decay_factor * phase1 * phase2;
		std::complex<double> ddyn_mult_dksi = ddecay_factor * phase1 * phase2 + decay_factor * phase1 * dphase2;
		
		double M_Q = 1.0 / one_minus_ksi;
		double M_C = -ksi / one_minus_ksi;
		
		Eigen::MatrixXcd static_grad = grad;
		
		grad(KSI, 0) = ddyn_mult_dksi * M_Q * (1.0 - eta) / 2.0 + dyn_mult * static_grad(KSI, 0);
		grad(KSI, 1) = ddyn_mult_dksi * M_C * (1.0 - eta) / 2.0 + dyn_mult * static_grad(KSI, 1);
		grad(KSI, 2) = ddyn_mult_dksi * M_C * (1.0 + eta) / 2.0 + dyn_mult * static_grad(KSI, 2);
		grad(KSI, 3) = ddyn_mult_dksi * M_Q * (1.0 + eta) / 2.0 + dyn_mult * static_grad(KSI, 3);
		
		grad(ETA, 0) = dyn_mult * static_grad(ETA, 0);
		grad(ETA, 1) = dyn_mult * static_grad(ETA, 1);
		grad(ETA, 2) = dyn_mult * static_grad(ETA, 2);
		grad(ETA, 3) = dyn_mult * static_grad(ETA, 3);
	}
	
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

