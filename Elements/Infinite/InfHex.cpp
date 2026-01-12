#include "InfHex.h"

double InfHex::compute_A() const {
	double cx = (x[0] + x[3] + x[4] + x[7]) / 4.0;
	double cy = (y[0] + y[3] + y[4] + y[7]) / 4.0;
	double cz = (z[0] + z[3] + z[4] + z[7]) / 4.0;
	
	double dx = cx - pole_x;
	double dy = cy - pole_y;
	double dz_val = cz - pole_z;
	
	return std::sqrt(dx*dx + dy*dy + dz_val*dz_val);
}

std::vector<std::complex<double>> InfHex::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> N;
	N.resize(8);
	
	double one_minus_ksi = 1.0 - ksi;
	double M0 = 1.0 / one_minus_ksi;
	double M1 = -ksi / one_minus_ksi;
	
	double N_eta_minus = (1.0 - eta) / 2.0;
	double N_eta_plus = (1.0 + eta) / 2.0;
	double N_zeta_minus = (1.0 - zeta) / 2.0;
	double N_zeta_plus = (1.0 + zeta) / 2.0;
	
	N[0] = M0 * N_eta_minus * N_zeta_minus;
	N[3] = M0 * N_eta_plus * N_zeta_minus;
	N[4] = M0 * N_eta_minus * N_zeta_plus;
	N[7] = M0 * N_eta_plus * N_zeta_plus;
	
	N[1] = M1 * N_eta_minus * N_zeta_minus;
	N[2] = M1 * N_eta_plus * N_zeta_minus;
	N[5] = M1 * N_eta_minus * N_zeta_plus;
	N[6] = M1 * N_eta_plus * N_zeta_plus;
	
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

Eigen::MatrixXcd InfHex::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd grad = Eigen::MatrixXcd::Zero(3, 8);
	
	double one_minus_ksi = 1.0 - ksi;
	double d = one_minus_ksi * one_minus_ksi;
	
	double dM0_dksi = 1.0 / d;
	double dM1_dksi = -(1.0 + ksi) / d;
	

	grad(KSI, 0) =  dM0_dksi * (1 - eta) * (1 - zeta) / 4.0;
	grad(KSI, 1) =  dM1_dksi * (1 - eta) * (1 - zeta) / 4.0;
	grad(KSI, 2) =  dM1_dksi * (1 + eta) * (1 - zeta) / 4.0;
	grad(KSI, 3) =  dM0_dksi * (1 + eta) * (1 - zeta) / 4.0;
	grad(KSI, 4) =  dM0_dksi * (1 - eta) * (1 + zeta) / 4.0;
	grad(KSI, 5) =  dM1_dksi * (1 - eta) * (1 + zeta) / 4.0;
	grad(KSI, 6) =  dM1_dksi * (1 + eta) * (1 + zeta) / 4.0;
	grad(KSI, 7) =  dM0_dksi * (1 + eta) * (1 + zeta) / 4.0;
	
	grad(ETA, 0) = -(1 - zeta) / (4 * one_minus_ksi);
	grad(ETA, 1) =  ksi * (1 - zeta) / (4 * one_minus_ksi);
	grad(ETA, 2) = -ksi * (1 - zeta) / (4 * one_minus_ksi);
	grad(ETA, 3) =  (1 - zeta) / (4 * one_minus_ksi);
	grad(ETA, 4) = -(1 + zeta) / (4 * one_minus_ksi);
	grad(ETA, 5) =  ksi * (1 + zeta) / (4 * one_minus_ksi);
	grad(ETA, 6) = -ksi * (1 + zeta) / (4 * one_minus_ksi);
	grad(ETA, 7) =  (1 + zeta) / (4 * one_minus_ksi);
	
	grad(ZETA, 0) = -(1 - eta) / (4 * one_minus_ksi);
	grad(ZETA, 1) =  ksi * (1 - eta) / (4 * one_minus_ksi);
	grad(ZETA, 2) =  ksi * (1 + eta) / (4 * one_minus_ksi);
	grad(ZETA, 3) = -(1 + eta) / (4 * one_minus_ksi);
	grad(ZETA, 4) =  (1 - eta) / (4 * one_minus_ksi);
	grad(ZETA, 5) = -ksi * (1 - eta) / (4 * one_minus_ksi);
	grad(ZETA, 6) = -ksi * (1 + eta) / (4 * one_minus_ksi);
	grad(ZETA, 7) =  (1 + eta) / (4 * one_minus_ksi);
	

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

		double M0 = 1.0 / one_minus_ksi;
		double M1 = -ksi / one_minus_ksi;
		
		double N_eta_minus = (1.0 - eta) / 2.0;
		double N_eta_plus = (1.0 + eta) / 2.0;
		double N_zeta_minus = (1.0 - zeta) / 2.0;
		double N_zeta_plus = (1.0 + zeta) / 2.0;
		
		Eigen::MatrixXcd static_grad = grad;
		
		grad(KSI, 0) = ddyn_mult_dksi * M0 * N_eta_minus * N_zeta_minus + dyn_mult * static_grad(KSI, 0);
		grad(KSI, 1) = ddyn_mult_dksi * M1 * N_eta_minus * N_zeta_minus + dyn_mult * static_grad(KSI, 1);
		grad(KSI, 2) = ddyn_mult_dksi * M1 * N_eta_plus * N_zeta_minus + dyn_mult * static_grad(KSI, 2);
		grad(KSI, 3) = ddyn_mult_dksi * M0 * N_eta_plus * N_zeta_minus + dyn_mult * static_grad(KSI, 3);
		grad(KSI, 4) = ddyn_mult_dksi * M0 * N_eta_minus * N_zeta_plus + dyn_mult * static_grad(KSI, 4);
		grad(KSI, 5) = ddyn_mult_dksi * M1 * N_eta_minus * N_zeta_plus + dyn_mult * static_grad(KSI, 5);
		grad(KSI, 6) = ddyn_mult_dksi * M1 * N_eta_plus * N_zeta_plus + dyn_mult * static_grad(KSI, 6);
		grad(KSI, 7) = ddyn_mult_dksi * M0 * N_eta_plus * N_zeta_plus + dyn_mult * static_grad(KSI, 7);
		
		for (int i = 0; i < 8; ++i) {
			grad(ETA, i) = dyn_mult * static_grad(ETA, i);
			grad(ZETA, i) = dyn_mult * static_grad(ZETA, i);
		}
	}
	
	return grad;
}

double InfHex::gaussPoint(LocVar var, int i) {
	static const double gr_0 = -1.0;
	static const double gr_1 = 1.0 / 3.0;
	
	static const double gp_std = 1.0 / std::sqrt(3.0);
	
	int zeta_idx = i / 4;
	int eta_idx = (i % 4) / 2;
	int ksi_idx = i % 2;
	
	if (var == KSI) {
		return (ksi_idx == 0) ? gr_0 : gr_1;
	}
	else if (var == ETA) {
		return (eta_idx == 0) ? -gp_std : gp_std;
	}
	else {
		return (zeta_idx == 0) ? -gp_std : gp_std;
	}
}

double InfHex::weight(LocVar var, int i) {
	
	int ksi_idx = i % 2;
	
	if (var == KSI) {  
		return (ksi_idx == 0) ? 0.5 : 1.5;
	}
	else {
		return 1.0;
	}
}

