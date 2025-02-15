#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

class logger {
	std::string file_name;
	std::ofstream file;

	logger() : file_name("log.txt"), file(file_name) {} // error open file
	~logger() { file.close(); }

	logger(const logger&) = delete;
	logger& operator=(const logger&) = delete;

public:
	static logger& log() {
		static logger Log;
		return Log;
	}
	void print(std::string msg) {
		std::time_t t = std::time(nullptr);
		struct tm timeInfo;
		errno_t err = localtime_s(&timeInfo, &t);
		char buffer[80];
		std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S", &timeInfo);

		file << buffer << "\t";
		file << msg << std::endl;
	}
};

