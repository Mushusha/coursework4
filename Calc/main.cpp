#include "Parser.h"
#include "Data.h"

#include "Tests.h"
#include "log1.h"

int main(int argc, char* argv[]) {
	logger &log = logger::log();
	log.print("Start program");
	//std::string file = argv[1];//"kqi400.fc";
	//std::shared_ptr<Parser> p = std::make_shared<Parser>();
	//p->read(file);
	//Data data(p);
	//data.set_output_param(std::vector<double>{0.1, 0.1}, std::vector<double>{5, 5}, 50);
	//data.solve();
	
	KirshError(argv[1], 1000000, 0.25);
	log.print("Done");
	
	//tests();
	return 0;
}


