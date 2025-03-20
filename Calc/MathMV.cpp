#include "MathMV.h"

void printMatrix(Eigen::SparseMatrix<double> A) {
	for (int i = 0; i < A.innerSize(); i++)
		for (int j = 0; j < A.outerSize(); j++)
			if (A.coeffRef(i, j) != 0.0)
				std::cout << i << " " << j << "\t" << A.coeffRef(i, j) << " " << std::endl;
}

void printVector(Eigen::SparseVector<double> B) {
	for (int i = 0; i < B.size(); i++)
		if (B.coeffRef(i) != 0.0)
			std::cout << i << "\t" << B.coeffRef(i) << std::endl;
}

double line(std::vector<double> A, std::vector<double> B, std::vector<double> X) {
	// return value y - f(x)
	return (X[1] - A[1]) * (B[0] - A[0]) - (X[0] - A[0]) * (B[1] - A[1]);
}

void output_points(std::vector<double> start, std::vector<double> end, int count, std::vector<std::vector<double>>& points) {
	points.resize(count);
	for (int i = 0; i < count; i++) {
		points[i].resize(start.size());
		for (int j = 0; j < start.size(); j++)
			points[i][j] = (end[j] - start[j]) * i / (count - 1) + start[j];
	}
}

//void KirshError(std::string filename, double P, double R) {
//	logger& log = logger::log();
//	log.print("Start calculate Kirsh error");
//
//	std::string file = filename;
//	std::shared_ptr<Parser> p = std::make_shared<Parser>();
//	p->read(file);
//	Data data(p);
//
//	data.set_output_param(std::vector<double>{0.8 * R, 0.8 * R}, std::vector<double>{5, 5}, 200);
//	//data.solve();
//
//	std::string fieldname;
//
//	std::vector<std::vector<double>> points;
//	output_points(data.line_start, data.line_end, data.points_count, points);
//	std::vector<double> points_x;
//	points_x.resize(data.points_count);
//	for (int i = 0; i < data.points_count; i++)
//		points_x[i] = points[i][0];
//
//
//	for (int comp = 0; comp < 3; comp++) {
//		field_name(fieldname, STRESS, comp, 2);
//
//		std::vector<double> analitics;
//		analitics.resize(points_x.size());
//
//		double an_max;
//		switch (comp) {
//		case Comp2D::XX:
//			for (int i = 0; i < points_x.size(); i++)
//				analitics[i] = P * (3 * pow(R / points_x[i], 4) / 4 - pow(R / points_x[i], 2)) / 2;
//			an_max = P / 6;
//			break;
//		case Comp2D::YY:
//			for (int i = 0; i < points_x.size(); i++)
//				analitics[i] = P * (-3 * pow(R / points_x[i], 4) / 4 + pow(R / points_x[i], 2) + 2) / 2;
//			an_max = P / 6 + P;
//			break;
//		case Comp2D::XY:
//			for (int i = 0; i < points_x.size(); i++)
//				analitics[i] = -P * pow(R / points_x[i], 2) / 4;
//			an_max = -P * pow(R / points_x[0], 2) / 4;
//			break;
//		default:
//			return;
//		}
//
//		std::vector<double> error;
//		error.resize(points_x.size());
//		for (int i = 0; i < points_x.size(); i++)
//			error[i] = data.out_stress[comp][i] - analitics[i];
//
//		auto compareAbs = [](const double& a, const double& b) {
//			return std::abs(a) < std::abs(b);
//		};
//		std::cout << fieldname + " error: " << std::abs(*std::max_element(error.begin(), error.end(), compareAbs) / an_max) * 100 << "%" << std::endl;
//		log.print("Kirsh error " + fieldname);
//	}
//}

//void meshConvergence(std::string name1, std::string name2, std::string name3) {
//	logger& log = logger::log();
//	log.print("Start calculate mesh convergence");
//
//	int count = 150;
//
//	std::string file1 = name1;
//	std::shared_ptr<Parser> p1 = std::make_shared<Parser>();
//	p1->read(file1);
//	Data data1(p1);
//	data1.set_output_param(std::vector<double>{0.3, -0.3}, std::vector<double>{3.9, -0.3}, count);
//	//data1.solve();
//
//	std::string file2 = name2;
//	std::shared_ptr<Parser> p2 = std::make_shared<Parser>();
//	p2->read(file2);
//	Data data2(p2);
//	data2.set_output_param(std::vector<double>{0.3, -0.3}, std::vector<double>{3.9, -0.3}, count);
//	//data2.solve();
//
//	std::string file3 = name3;
//	std::shared_ptr<Parser> p3 = std::make_shared<Parser>();
//	p3->read(file3);
//	Data data3(p3);
//	data3.set_output_param(std::vector<double>{0.1, -0.3}, std::vector<double>{3.9, -0.3}, count);
//	//data3.solve();
//
//	std::string fieldname;
//
//	auto compareAbs = [](const double& a, const double& b) {
//		return std::abs(a) < std::abs(b);
//	};
//
//	for (int comp = 0; comp < 3; comp++) {
//		field_name(fieldname, STRESS, comp, 2);
//
//		std::vector<double> error1;
//		error1.resize(count);
//		for (int i = 0; i < count; i++)
//			error1[i] = data1.out_stress[comp][i] - data2.out_stress[comp][i];
//		std::vector<double> error2;
//		error2.resize(count);
//		for (int i = 0; i < count; i++)
//			error2[i] = data2.out_stress[comp][i] - data3.out_stress[comp][i];
//
//		std::cout << std::endl;
//		std::cout << fieldname + " error1: " << std::abs(*std::max_element(error1.begin(), error1.end(), compareAbs) / *std::max_element(data1.out_stress[comp].begin(), data1.out_stress[comp].end(), compareAbs)) * 100 << "%" << std::endl;
//		std::cout << fieldname + " error2: " << std::abs(*std::max_element(error2.begin(), error2.end(), compareAbs) / *std::max_element(data2.out_stress[comp].begin(), data2.out_stress[comp].end(), compareAbs)) * 100 << "%" << std::endl;
//		log.print("mesh convergence " + fieldname);
//	}
//}

