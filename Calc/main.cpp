#include "Calc.h"

#include "Tests/Tests.h"
#include "log1.h"

int main(int argc, char* argv[]) {
	Test();

	logger &log = logger::log();
	log.print("Start program");

	if (argv[1]) {
		Calculate solve(argv[1]);
		solve.Solve();
	}

	log.print("Done");
	
	// tests();
	return 0;
}


