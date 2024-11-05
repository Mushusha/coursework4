#include "Parser.h"
#include "Data.h"

int main(int argc, char* argv[])
{
	// time
	string file = "C:/Users/mushu/Desktop/cw/forFC/k.fc";
	std::shared_ptr<Parser> p = std::make_shared<Parser>();
	p->read(file);
	Data data(p);
	std::cout << 1;
	//Eigen::SparseMatrix <double> K(p->mesh.elems_count * 2, p->mesh.elems_count * 2);
	//K = globalK(p, 4, 2);
	//for (int k = 0; k < K.outerSize() - 1000; ++k)
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it)
	//	{
	//		std::cout << it.value() << "\t---" << it.row() << "\t--" << it.col() << endl;
	//	}
	return 0;
}
