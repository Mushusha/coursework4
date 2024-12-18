#include "Parser.h"
#include "Data.h"

#include "Tests.h"
#include "log.h"

int main(int argc, char* argv[]) {
	logger &log = logger::log();
	log.print("Start program");
	//std::string file = "k.fc";
	//std::shared_ptr<Parser> p = std::make_shared<Parser>();
	//p->read(file);
	//Data data(p);
	//data.solve();
	////data.set_output_param(std::vector<double>{0.0, -4.5}, std::vector<double>{0.0, 10}, 100, 1000000, 0.25); 
	//data.set_output_param(std::vector<double>{0.85, 0.85}, std::vector<double>{50, 50}, 300);
	//data.smoothing();

	//KirshError("k21468.fc", 1, 1);
	meshConvergence("kot5835.fc", "kot11614.fc", "kot23606.fc");
	log.print("Done");

	//tests();
	return 0;
}
