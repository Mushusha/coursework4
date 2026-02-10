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
#include "Data.h"

template <int NODES>
class SpectralHex : public Element {
	static_assert(NODES >= 2, "NODES must be >= 2");

public:
	SpectralHex() : Element() {
		init_gll();
	}
	SpectralHex(int id, ElemType type, std::vector<int> nodes)
		: Element(id, type, nodes) {
		init_gll();
	}
	SpectralHex(const SpectralHex& other) : Element(other) {
		init_gll();
	}
	SpectralHex& operator=(const SpectralHex& other) {
		if (this != &other) Element::operator=(other);
		return *this;
	}
	SpectralHex(SpectralHex&& other) noexcept : Element(std::move(other)) {
		init_gll();
	}
	SpectralHex& operator=(SpectralHex&& other) noexcept {
		if (this != &other) Element::operator=(std::move(other));
		return *this;
	}
	virtual ~SpectralHex() = default;

	std::vector<std::complex<double>> FF(double ksi, double eta, double zeta) override;

	Eigen::MatrixXcd B(double ksi = 0, double eta = 0, double zeta = 0) override;

	Eigen::MatrixXcd localK() override;
	std::vector<double> localF(double mult = 1) override;

	Eigen::MatrixXd localC() override;
	std::vector<double> localR(std::vector<double> value) override;
	Eigen::MatrixXcd localM() override;
	Eigen::MatrixXd localDamping() override;

	std::vector<int> edge_to_node(int edge) override;

	double gaussPoint(LocVar var, int i) override;
	double weight(LocVar var, int i) override;
	double Volume() final;

	Eigen::MatrixXcd gradFF(double ksi, double eta, double zeta) override;
	Eigen::MatrixXcd J(double ksi, double eta, double zeta) override;

	bool pointInElem(std::vector<double> point) override;
	std::vector<double> coordFF(double x0, double y0, double z0) override;

protected:
	void set_pressure(int edge, double value) override;
	
	std::vector<double> gll_x;
	std::vector<double> gll_w;

	double lagrange1D(int i, double x) const;
	double dlagrange1D(int i, double x) const;

	void init_gll();

private:
	static constexpr int nNodes = NODES * NODES * NODES;

	double area_face(int face);
	std::array<double, 3> normal(int face);

	void get_face_center(int face, double& ksi, double& eta, double& zeta);
};

template<int NODES>
void SpectralHex<NODES>::init_gll() {
	if (gll_x.empty() || gll_w.empty())
		compute_gll_nodes_weights(NODES - 1, gll_x, gll_w);
}

template<int NODES>
double SpectralHex<NODES>::lagrange1D(int i, double x) const {
	double val = 1.0;
	for (int m = 0; m < NODES; ++m) {
		if (m == i) continue;
		val *= (x - gll_x[m]) / (gll_x[i] - gll_x[m]);
	}
	return val;
}

template<int NODES>
double SpectralHex<NODES>::dlagrange1D(int i, double x) const {
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
std::vector<std::complex<double>> SpectralHex<NODES>::FF(double ksi, double eta, double zeta) {
	std::vector<std::complex<double>> N(nNodes);
	int idx = 0;

	for (int k = 0; k < NODES; ++k) {
		double Nk_zeta = lagrange1D(k, zeta);
		for (int j = 0; j < NODES; ++j) {
			double Nj_eta = lagrange1D(j, eta);
			for (int i = 0; i < NODES; ++i) {
				double Ni_ksi = lagrange1D(i, ksi);
				double value = Ni_ksi * Nj_eta * Nk_zeta;
				N[idx] = std::complex<double>(value, 0.0);
				idx++;
			}
		}
	}
	return N;
}

template<int NODES>
Eigen::MatrixXcd SpectralHex<NODES>::gradFF(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd grad(3, nNodes);
	grad.setZero();
	int idx = 0;

	for (int k = 0; k < NODES; ++k) {
		double Nk_zeta = lagrange1D(k, zeta);
		double dNk_dzeta = dlagrange1D(k, zeta);

		for (int j = 0; j < NODES; ++j) {
			double Nj_eta = lagrange1D(j, eta);
			double dNj_deta = dlagrange1D(j, eta);

			for (int i = 0; i < NODES; ++i, ++idx) {
				double Ni_ksi = lagrange1D(i, ksi);
				double dNi_dksi = dlagrange1D(i, ksi);

				grad(0, idx) = std::complex<double>(dNi_dksi * Nj_eta * Nk_zeta, 0.0);
				grad(1, idx) = std::complex<double>(Ni_ksi * dNj_deta * Nk_zeta, 0.0);
				grad(2, idx) = std::complex<double>(Ni_ksi * Nj_eta * dNk_dzeta, 0.0);
			}
		}
	}
	return grad;
}

template<int NODES>
Eigen::MatrixXcd SpectralHex<NODES>::J(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd jac = Eigen::MatrixXcd::Zero(3, 3);
	Eigen::MatrixXcd g = gradFF(ksi, eta, zeta);

	for (int a = 0; a < nNodes; ++a) {
		jac(0, 0) += g(KSI, a) * x[a];
		jac(0, 1) += g(KSI, a) * y[a];
		jac(0, 2) += g(KSI, a) * z[a];
		jac(1, 0) += g(ETA, a) * x[a];
		jac(1, 1) += g(ETA, a) * y[a];
		jac(1, 2) += g(ETA, a) * z[a];
		jac(2, 0) += g(ZETA, a) * x[a];
		jac(2, 1) += g(ZETA, a) * y[a];
		jac(2, 2) += g(ZETA, a) * z[a];
	}
	return jac;
}

template<int NODES>
Eigen::MatrixXcd SpectralHex<NODES>::B(double ksi, double eta, double zeta) {
	Eigen::MatrixXcd Bm = Eigen::MatrixXcd::Zero(6, 3 * nNodes);
	Eigen::Matrix3cd invJ = J(ksi, eta, zeta).inverse();
	Eigen::MatrixXcd dN = invJ * gradFF(ksi, eta, zeta);

	for (int a = 0; a < nNodes; ++a) {
		double dNx = dN(X, a).real();
		double dNy = dN(Y, a).real();
		double dNz = dN(Z, a).real();

		Bm(0, 3 * a) = dNx;
		Bm(1, 3 * a + 1) = dNy;
		Bm(2, 3 * a + 2) = dNz;

		Bm(3, 3 * a) = dNy;
		Bm(3, 3 * a + 1) = dNx;

		Bm(4, 3 * a + 1) = dNz;
		Bm(4, 3 * a + 2) = dNy;

		Bm(5, 3 * a) = dNz;
		Bm(5, 3 * a + 2) = dNx;
	}
	return Bm;
}

template<int NODES>
Eigen::MatrixXcd SpectralHex<NODES>::localK() {
	Eigen::MatrixXcd K = Eigen::MatrixXcd::Zero(3 * nNodes, 3 * nNodes);

	for (int k = 0; k < NODES; ++k) {
		for (int j = 0; j < NODES; ++j) {
			for (int i = 0; i < NODES; ++i) {
				int node_idx = i + j * NODES + k * NODES * NODES;
				double ksi = gaussPoint(KSI, node_idx);
				double eta = gaussPoint(ETA, node_idx);
				double zeta = gaussPoint(ZETA, node_idx);
				double w_ksi = weight(KSI, node_idx);
				double w_eta = weight(ETA, node_idx);
				double w_zeta = weight(ZETA, node_idx);

				Eigen::MatrixXcd Bm = B(ksi, eta, zeta);
				Eigen::MatrixXcd Jm = J(ksi, eta, zeta);

				double detJ = std::abs(Jm.determinant());

				K += (w_ksi * w_eta * w_zeta) * (Bm.adjoint() * D.template cast<std::complex<double>>() * Bm * detJ);
			}
		}
	}
	return K;
}

template<int NODES>
std::vector<double> SpectralHex<NODES>::localF(double mult) {
	std::vector<double> F(3 * nNodes, 0.0);
	if (load.empty()) return F;

	for (auto const& l : load) {
		const int face = l.first.first;
		int comp = l.first.second;
		double value = l.second;

		std::vector<int> face_nodes = edge_to_node(face);
		int m = face_nodes.size();

		for (int i = 0; i < NODES; ++i) {
			double u = gll_x[i];
			double w_u = gll_w[i];

			for (int j = 0; j < NODES; ++j) {
				double v = gll_x[j];
				double w_v = gll_w[j];

				double ksi, eta, zeta;
				// Порядок как в сетке: 0=zeta=-1, 1=eta=-1, 2=ksi=+1, 3=eta=+1, 4=ksi=-1, 5=zeta=+1
				switch (face) {
				case 0: ksi = u; eta = v; zeta = -1.0; break;
				case 1: ksi = u; eta = -1.0; zeta = v; break;
				case 2: ksi = 1.0; eta = u; zeta = v; break;
				case 3: ksi = u; eta = 1.0; zeta = v; break;
				case 4: ksi = -1.0; eta = u; zeta = v; break;
				case 5: ksi = u; eta = v; zeta = 1.0; break;
				default: ksi = u; eta = v; zeta = -1.0; break;
				}

				Eigen::MatrixXcd grad = gradFF(ksi, eta, zeta);
				Eigen::Vector3d dX_du(0, 0, 0), dX_dv(0, 0, 0);
				// 0,5 — zeta; 1,3 — eta; 2,4 — ksi
				switch (face) {
				case 0: case 5:
					for (int a = 0; a < nNodes; ++a) {
						dX_du[0] += grad(0, a).real() * x[a]; dX_du[1] += grad(0, a).real() * y[a]; dX_du[2] += grad(0, a).real() * z[a];
						dX_dv[0] += grad(1, a).real() * x[a]; dX_dv[1] += grad(1, a).real() * y[a]; dX_dv[2] += grad(1, a).real() * z[a];
					}
					break;
				case 1: case 3:
					for (int a = 0; a < nNodes; ++a) {
						dX_du[0] += grad(0, a).real() * x[a]; dX_du[1] += grad(0, a).real() * y[a]; dX_du[2] += grad(0, a).real() * z[a];
						dX_dv[0] += grad(2, a).real() * x[a]; dX_dv[1] += grad(2, a).real() * y[a]; dX_dv[2] += grad(2, a).real() * z[a];
					}
					break;
				case 2: case 4:
					for (int a = 0; a < nNodes; ++a) {
						dX_du[0] += grad(1, a).real() * x[a]; dX_du[1] += grad(1, a).real() * y[a]; dX_du[2] += grad(1, a).real() * z[a];
						dX_dv[0] += grad(2, a).real() * x[a]; dX_dv[1] += grad(2, a).real() * y[a]; dX_dv[2] += grad(2, a).real() * z[a];
					}
					break;
				default: break;
				}

				Eigen::Vector3d normal_vec = dX_du.cross(dX_dv);
				double dS = normal_vec.norm() * w_u * w_v;

				auto Nvals = FF(ksi, eta, zeta);

				for (int local_node_index : face_nodes) {
					double contribution = mult * value * Nvals[local_node_index].real() * dS;
					F[3 * local_node_index + comp] += contribution;
				}
			}
		}
	}

	return F;
}

template<int NODES>
Eigen::MatrixXd SpectralHex<NODES>::localC() {
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(nNodes, nNodes);

	for (int k = 0; k < NODES; ++k) {
		for (int j = 0; j < NODES; ++j) {
			for (int i = 0; i < NODES; ++i) {
				int node_idx = i + j * NODES + k * NODES * NODES;
				double ksi = gaussPoint(KSI, node_idx);
				double eta = gaussPoint(ETA, node_idx);
				double zeta = gaussPoint(ZETA, node_idx);
				double w_ksi = weight(KSI, node_idx);
				double w_eta = weight(ETA, node_idx);
				double w_zeta = weight(ZETA, node_idx);

				double detJ = std::abs(J(ksi, eta, zeta).determinant());
				auto Nvals = FF(ksi, eta, zeta);

				for (int a = 0; a < nNodes; ++a)
					for (int b = 0; b < nNodes; ++b)
						C(a, b) += w_ksi * w_eta * w_zeta * Nvals[a].real() * Nvals[b].real() * detJ;
			}
		}
	}
	return C;
}

template<int NODES>
Eigen::MatrixXd SpectralHex<NODES>::localDamping() {
	if (density <= 0.0) 
		return Eigen::MatrixXd::Zero(nNodes, nNodes);
	return Data::damping * density * localC();
}

template<int NODES>
std::vector<double> SpectralHex<NODES>::localR(std::vector<double> value) {
	std::vector<double> R(nNodes, 0.0);

	for (int k = 0; k < NODES; ++k) {
		for (int j = 0; j < NODES; ++j) {
			for (int i = 0; i < NODES; ++i) {
				int node_idx = i + j * NODES + k * NODES * NODES;
				double ksi = gaussPoint(KSI, node_idx);
				double eta = gaussPoint(ETA, node_idx);
				double zeta = gaussPoint(ZETA, node_idx);
				double w_ksi = weight(KSI, node_idx);
				double w_eta = weight(ETA, node_idx);
				double w_zeta = weight(ZETA, node_idx);

				double detJ = std::abs(J(ksi, eta, zeta).determinant());
				auto Nvals = FF(ksi, eta, zeta);

				for (int a = 0; a < nNodes; ++a)
					R[a] += value[a] * w_ksi * w_eta * w_zeta * Nvals[a].real() * detJ;
			}
		}
	}
	return R;
}

template<int NODES>
Eigen::MatrixXcd SpectralHex<NODES>::localM() {
	if (density == 0.0)
		throw std::runtime_error("Error: density is zero in element " + std::to_string(id));

	Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(3 * nNodes, 3 * nNodes);

	for (int k = 0; k < NODES; ++k) {
		for (int j = 0; j < NODES; ++j) {
			for (int i = 0; i < NODES; ++i) {
				int node_idx = i + j * NODES + k * NODES * NODES;
				double ksi = gaussPoint(KSI, node_idx);
				double eta = gaussPoint(ETA, node_idx);
				double zeta = gaussPoint(ZETA, node_idx);
				double w_ksi = weight(KSI, node_idx);
				double w_eta = weight(ETA, node_idx);
				double w_zeta = weight(ZETA, node_idx);

				double detJ = std::abs(J(ksi, eta, zeta).determinant());
				auto Nvals = FF(ksi, eta, zeta);

				for (int a = 0; a < nNodes; ++a)
					for (int b = 0; b < nNodes; ++b) {
						std::complex<double> Mij = w_ksi * w_eta * w_zeta * std::conj(Nvals[a]) * Nvals[b] * detJ;
						M(3 * a, 3 * b) += Mij;
						M(3 * a + 1, 3 * b + 1) += Mij;
						M(3 * a + 2, 3 * b + 2) += Mij;
					}
			}
		}
	}
	return density * M;
}

template<int NODES>
std::vector<int> SpectralHex<NODES>::edge_to_node(int face) {
	std::vector<int> res;
	// Порядок граней как в сетке (Gmsh/solver): 0=zeta=-1, 1=eta=-1, 2=ksi=+1, 3=eta=+1, 4=ksi=-1, 5=zeta=+1
	// idx = i + j*NODES + k*NODES*NODES (i=ksi, j=eta, k=zeta)
	switch (face) {
	case 0: // mesh 0 = zeta = -1, k = 0
		for (int i = 0; i < NODES; ++i) {
			for (int j = 0; j < NODES; ++j) {
				res.push_back(i + j * NODES + 0 * NODES * NODES);
			}
		}
		break;
	case 1: // mesh 1 = eta = -1, j = 0
		for (int i = 0; i < NODES; ++i) {
			for (int k = 0; k < NODES; ++k) {
				res.push_back(i + 0 * NODES + k * NODES * NODES);
			}
		}
		break;
	case 2: // mesh 2 = ksi = +1, i = NODES-1
		for (int k = 0; k < NODES; ++k) {
			for (int j = 0; j < NODES; ++j) {
				res.push_back((NODES - 1) + j * NODES + k * NODES * NODES);
			}
		}
		break;
	case 3: // mesh 3 = eta = +1, j = NODES-1
		for (int k = 0; k < NODES; ++k) {
			for (int i = 0; i < NODES; ++i) {
				res.push_back(i + (NODES - 1) * NODES + k * NODES * NODES);
			}
		}
		break;
	case 4: // mesh 4 = ksi = -1, i = 0
		for (int k = 0; k < NODES; ++k) {
			for (int j = 0; j < NODES; ++j) {
				res.push_back(0 + j * NODES + k * NODES * NODES);
			}
		}
		break;
	case 5: // mesh 5 = zeta = +1, k = NODES-1
		for (int i = 0; i < NODES; ++i) {
			for (int j = 0; j < NODES; ++j) {
				res.push_back(i + j * NODES + (NODES - 1) * NODES * NODES);
			}
		}
		break;
	default:
		throw std::runtime_error("Wrong face index in edge_to_node");
	}

	return res;
}

template<int NODES>
double SpectralHex<NODES>::gaussPoint(LocVar var, int i) {
	if (i < 0 || i >= NODES * NODES * NODES)
		throw std::out_of_range("SpectralHex<NODES>::gaussPoint: index out of range");

	int idx_zeta = i / (NODES * NODES);
	int idx_eta = (i % (NODES * NODES)) / NODES;
	int idx_ksi = i % NODES;

	if (var == KSI)
		return gll_x[idx_ksi];
	else if (var == ETA)
		return gll_x[idx_eta];
	else if (var == ZETA)
		return gll_x[idx_zeta];
	else
		return 0.0;
}

template<int NODES>
double SpectralHex<NODES>::weight(LocVar var, int i) {
	if (i < 0 || i >= NODES * NODES * NODES)
		throw std::out_of_range("SpectralHex<NODES>::weight: index out of range");

	int idx_zeta = i / (NODES * NODES);
	int idx_eta = (i % (NODES * NODES)) / NODES;
	int idx_ksi = i % NODES;

	if (var == KSI)
		return gll_w[idx_ksi];
	else if (var == ETA)
		return gll_w[idx_eta];
	else if (var == ZETA)
		return gll_w[idx_zeta];
	else
		return 0.0;
}

template<int NODES>
double SpectralHex<NODES>::Volume() {
	std::complex<double> V = 0.0;

	for (int k = 0; k < NODES; ++k) {
		double zeta = gll_x[k];
		double w_zeta = gll_w[k];

		for (int j = 0; j < NODES; ++j) {
			double eta = gll_x[j];
			double w_eta = gll_w[j];

			for (int i = 0; i < NODES; ++i) {
				double ksi = gll_x[i];
				double w_ksi = gll_w[i];

				auto Nvals = FF(ksi, eta, zeta);
				double detJ = std::abs(J(ksi, eta, zeta).determinant());

				for (int a = 0; a < nNodes; ++a)
					V += w_ksi * w_eta * w_zeta * Nvals[a] * detJ;
			}
		}
	}
	return V.real();
}

template<int NODES>
bool SpectralHex<NODES>::pointInElem(std::vector<double> point) {
	return false;
}

template<int NODES>
std::vector<double> SpectralHex<NODES>::coordFF(double x0, double y0, double z0) {
	return std::vector<double>();
}

template<int NODES>
void SpectralHex<NODES>::set_pressure(int face, double value) {
	std::array<double, 3> n = normal(face);

	for (int i = 0; i < 3; i++) {
		std::pair <int, int> pair(face, i);
		load.insert({ pair, value * n[i] });
	}

}

template<int NODES>
double SpectralHex<NODES>::area_face(int face) {
	double area = 0.0;
	int m = NODES;

	for (int i = 0; i < m; ++i) {
		double u = gll_x[i];
		double w_u = gll_w[i];

		for (int j = 0; j < m; ++j) {
			double v = gll_x[j];
			double w_v = gll_w[j];

			double ksi, eta, zeta;

			switch (face) {
			case 0: ksi = u; eta = v; zeta = -1.0; break;
			case 1: ksi = -1.0; eta = u; zeta = v; break;
			case 2: ksi = 1.0; eta = u; zeta = v; break;
			case 3: ksi = u; eta = 1.0; zeta = v; break;
			case 4: ksi = u; eta = u; zeta = 1.0; break;
			case 5: ksi = u; eta = -1.0; zeta = v; break;
			}

			Eigen::MatrixXcd grad = gradFF(ksi, eta, zeta);
			Eigen::Vector3d dX_du(0, 0, 0), dX_dv(0, 0, 0);

			switch (face) {
			case 5: case 3:
				for (int a = 0; a < nNodes; ++a) {
					dX_du[0] += grad(0, a).real() * x[a];
					dX_du[1] += grad(0, a).real() * y[a];
					dX_du[2] += grad(0, a).real() * z[a];

					dX_dv[0] += grad(2, a).real() * x[a];
					dX_dv[1] += grad(2, a).real() * y[a];
					dX_dv[2] += grad(2, a).real() * z[a];
				}
				break;
			case 4: case 0:
				for (int a = 0; a < nNodes; ++a) {
					dX_du[0] += grad(0, a).real() * x[a];
					dX_du[1] += grad(0, a).real() * y[a];
					dX_du[2] += grad(0, a).real() * z[a];

					dX_dv[0] += grad(1, a).real() * x[a];
					dX_dv[1] += grad(1, a).real() * y[a];
					dX_dv[2] += grad(1, a).real() * z[a];
				}
				break;
			case 2: case 1:
				for (int a = 0; a < nNodes; ++a) {
					dX_du[0] += grad(1, a).real() * x[a];
					dX_du[1] += grad(1, a).real() * y[a];
					dX_du[2] += grad(1, a).real() * z[a];

					dX_dv[0] += grad(2, a).real() * x[a];
					dX_dv[1] += grad(2, a).real() * y[a];
					dX_dv[2] += grad(2, a).real() * z[a];
				}
				break;
			}

			Eigen::Vector3d normal_vec = dX_du.cross(dX_dv);
			double dS = normal_vec.norm() * w_u * w_v;
			area += dS;
		}
	}

	return area;
}

template<int NODES>
std::array<double, 3> SpectralHex<NODES>::normal(int face) {
	double ksi, eta, zeta;
	get_face_center(face, ksi, eta, zeta);

	Eigen::MatrixXcd Jm = J(ksi, eta, zeta);
	Eigen::Vector3d dX_dksi = Jm.col(0).real();
	Eigen::Vector3d dX_deta = Jm.col(1).real();
	Eigen::Vector3d dX_dzeta = Jm.col(2).real();

	Eigen::Vector3d normal_vec;

	switch (face) {
	case 0:
		normal_vec = dX_dksi.cross(dX_dzeta);
		break;
	case 1:
		normal_vec = dX_dzeta.cross(dX_deta);
		break;
	case 2:
		normal_vec = dX_deta.cross(dX_dzeta);
		break;
	case 3:
		normal_vec = dX_dzeta.cross(dX_dksi);
		break;
	case 4:
		normal_vec = dX_dzeta.cross(dX_deta);
		break;
	case 5:
		normal_vec = dX_dksi.cross(dX_deta);
		break;
	default:
		throw std::runtime_error("Wrong face index in normal");
	}

	normal_vec.normalize();
	return { normal_vec[0], normal_vec[1], normal_vec[2] };
}

template<int NODES>
void SpectralHex<NODES>::get_face_center(int face, double& ksi, double& eta, double& zeta) {
	switch (face) {
	case 0: ksi = 0.0; eta = -1.0; zeta = 0.0; break;
	case 1: ksi = -1.0; eta = 0.0; zeta = 0.0; break;
	case 2: ksi = 1.0; eta = 0.0; zeta = 0.0; break;
	case 3: ksi = 0.0; eta = 1.0; zeta = 0.0; break;
	case 4: ksi = -1.0; eta = 0.0; zeta = 0.0; break;
	case 5: ksi = 0.0; eta = 0.0; zeta = -1.0; break;
	default:
		throw std::runtime_error("Wrong face index in get_face_center");
	}
}