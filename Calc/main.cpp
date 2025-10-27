#include <iostream>
#include <csignal>
#include <cstdlib>

#include "Server.h"
#include "Tests/Tests.h"
#include "log1.h"


std::unique_ptr<ApiServer> server;

void signalHandler(int signal) {
	std::exit(0);
}

int main(int argc, char* argv[]) {
	Test();

	logger& log = logger::log();

	std::signal(SIGINT, signalHandler);
	std::signal(SIGTERM, signalHandler);

	try {
		if (argv[1]) {
			logger& log = logger::log();
			log.print("Start program");

			log.print(argv[1]);
			Calculate solve(argv[1]);
			solve.Solve();

			log.print("Done");

			return 0;
		}

		server = std::make_unique<ApiServer>(3000);
		server->run();
	}
	catch (const std::exception& e) {
		std::cerr << "FATAL ERROR: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}