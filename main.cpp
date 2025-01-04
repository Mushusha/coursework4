#include "Parser.h"
#include "Data.h"

#include "Tests.h"
#include "log.h"

int main(int argc, char* argv[]) {
	logger &log = logger::log();
	log.print("Start program");
	//std::string file = "k21468.fc";
	//std::shared_ptr<Parser> p = std::make_shared<Parser>();
	//p->read(file);
	//Data data(p);
	//data.set_output_param(std::vector<double>{0.0, 0.0}, std::vector<double>{5, 5}, 300);
	//data.solve();

	KirshError("k21468.fc", 1, 1);

	log.print("Done");

	//tests();
	return 0;
}
