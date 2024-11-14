#include "Parser.h"
#include "Data.h"

#include "Tests.h"
#include "log.h"

int main(int argc, char* argv[]) {

	// errors

	logger &log = logger::log();
	log.print("Start program");
	string file = "C:/Users/mushu/Desktop/cwgit/coursework4/k.fc";
	std::shared_ptr<Parser> p = std::make_shared<Parser>();
	p->read(file);
	Data data(p);
	data.solve();
	//std::cout << data.U << std::endl;
	log.print("Done");
	tests();
	return 0;
}
