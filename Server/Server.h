#pragma once
#include <crow.h>
#include <crow/middlewares/cors.h>

#include <string>
#include <memory>
#include <iostream>
#include <chrono>
#include <vector>

#include "Calc.h"
#include "log1.h"


struct VtuData {
    std::vector<double> points;
    std::vector<double> displacements;
    std::vector<double> stresses;
    std::vector<double> strains;
    std::vector<int> connectivity;
    std::vector<int> elementTypes;

    int numberOfPoints = 0;
    int numberOfCells = 0;

    crow::json::wvalue toJson() const;
};

class Calculate;
class ApiServer {
public:
    ApiServer(int port = 3000);
    void run();

private:
    void setupRoutes();
    VtuData parseVtuFile(const std::string& filename);
    std::string getVtuFilename(const std::string& fcFilename);

    crow::App<crow::CORSHandler> app_;
    int port;
    std::string curr_filename;
};