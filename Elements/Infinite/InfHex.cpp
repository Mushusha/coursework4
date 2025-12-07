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
		std::complex<double> phase2 = std::exp(i * k_wave * A_len / one_minus_ksi);
		std::complex<double> mult = decay_factor * phase1 * phase2;
		
		std::transform(N.begin(), N.end(), N.begin(), 
			[&mult](std::complex<double> a) { return a * mult; });
	}
	
	return N;
}

Eigen::MatrixXcd InfHex::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd grad = Eigen::MatrixXcd::Zero(3, 8);
	
	double d = (1 - ksi) * (1 - ksi);
	
	grad(KSI, 0) =  (1 - eta) * (1 - zeta) / (4 * d);
	grad(KSI, 1) = -(1 - eta) * (1 - zeta) / (4 * d);
	grad(KSI, 2) = -(1 + eta) * (1 - zeta) / (4 * d);
	grad(KSI, 3) =  (1 + eta) * (1 - zeta) / (4 * d);
	grad(KSI, 4) =  (1 - eta) * (1 + zeta) / (4 * d);
	grad(KSI, 5) = -(1 - eta) * (1 + zeta) / (4 * d);
	grad(KSI, 6) = -(1 + eta) * (1 + zeta) / (4 * d);
	grad(KSI, 7) =  (1 + eta) * (1 + zeta) / (4 * d);
	
	grad(ETA, 0) = -(1 - zeta) / (4 * (1 - ksi));
	grad(ETA, 1) =  ksi * (1 - zeta) / (4 * (1 - ksi));
	grad(ETA, 2) = -ksi * (1 - zeta) / (4 * (1 - ksi));
	grad(ETA, 3) =  (1 - zeta) / (4 * (1 - ksi));
	grad(ETA, 4) = -(1 + zeta) / (4 * (1 - ksi));
	grad(ETA, 5) =  ksi * (1 + zeta) / (4 * (1 - ksi));
	grad(ETA, 6) = -ksi * (1 + zeta) / (4 * (1 - ksi));
	grad(ETA, 7) =  (1 + zeta) / (4 * (1 - ksi));
	
	grad(ZETA, 0) = -(1 - eta) / (4 * (1 - ksi));
	grad(ZETA, 1) =  ksi * (1 - eta) / (4 * (1 - ksi));
	grad(ZETA, 2) =  ksi * (1 + eta) / (4 * (1 - ksi));
	grad(ZETA, 3) = -(1 + eta) / (4 * (1 - ksi));
	grad(ZETA, 4) =  (1 - eta) / (4 * (1 - ksi));
	grad(ZETA, 5) = -ksi * (1 - eta) / (4 * (1 - ksi));
	grad(ZETA, 6) = -ksi * (1 + eta) / (4 * (1 - ksi));
	grad(ZETA, 7) =  (1 + eta) / (4 * (1 - ksi));
	
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

