//#pragma once
//
//#include "QuadElements.h"
//#include "Element.h"
//
//
//class infQuadElement : public  quadElement {
//public:
//	infQuadElement() : quadElement() {}
//	infQuadElement(int id, int type, std::vector<int> nodes)
//		: quadElement(id, type, nodes) {
//	}
//	virtual ~infQuadElement() = default;
//
//	//Eigen::MatrixXd localK();
//	//std::vector <double> localF() override;
//	//Eigen::MatrixXd B(double ksi = 0, double eta = 0, double zeta = 0) override;
//	//std::vector<double> FF(double ksi, double eta, double zeta = 0);
//
//	//bool pointInElem(std::vector<double> point) override;
//
//	//Eigen::MatrixXd localC() override;
//	//std::vector <double> localR(std::vector<double> value) override;
//
//	//std::vector<double> coordFF(double x0, double y0, double z0 = 0) override;
//
////protected:
////	Eigen::MatrixXd  gradFF(double ksi, double eta, double zeta = 0) override;
////	Eigen::MatrixXd J(double ksi, double eta, double zeta = 0) override;
////
////	void set_pressure(int edge, double value);
////private:
////	double len_edge(int edge);
////	std::pair<int, int> edge_to_node(int edge);
//};

