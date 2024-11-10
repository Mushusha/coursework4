#include "Parser.h"
#include "Data.h"

int main(int argc, char* argv[])
{
	// time
	// log
	// errors
	// tests

	string file = "C:/Users/mushu/Desktop/cwgit/coursework4/k.fc";
	std::shared_ptr<Parser> p = std::make_shared<Parser>();
	p->read(file);
	Data data(p);
	data.solve();
	std::cout << data.U << std::endl;
	std::cout << "done";

	return 0;
}
