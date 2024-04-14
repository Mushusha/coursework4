#include "Parser.h"
#include "global.h"
#include "QuadElements.h"

int main(int argc, char* argv[])
{
	string file = "C:/Users/mushu/Desktop/cw/forFC/cw.fc";
	std::shared_ptr<Parser> P = std::make_shared<Parser>();
	P->read(file);
	//for (int i = 0; i < P->nodesets.size(); i++)
	//	for (int j = 0; j < P->nodesets[i].apply_to.size(); j++)
	//		std::cout << P->nodesets[i].apply_to[j] << endl;
	////std::cout << P->mesh.elems_count;
	std::vector <std::shared_ptr <Element>> elem;
	elem = createElements(P, 4);
	return 0;
}


//int main() {
//
//	vector <float> x = { -1, 0, 1, 2, 1, 3, 7, -8 };
//	vector <float> y = { 2, 1, -2, 0, 0, -5, 6, 10 };
//	vector <float> z = { 0, -1, 1, 0, 6, 5, 3, -2 };
//
//	cout << hexaLocK(x, y, z, hexaInfN);
//
//	return 0;
//}