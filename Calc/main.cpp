#include "Calc.h"

#include "Tests.h"
#include "log1.h"

int main(int argc, char* argv[]) {
	logger &log = logger::log();
	log.print("Start program");

	Calculate solve(argv[1]);
	solve.Solve();

	log.print("Done");
	
	// tests();
	return 0;
}


