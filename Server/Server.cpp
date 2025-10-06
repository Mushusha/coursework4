#include "Server.h"

#include <fstream>
#include <sstream>
#include <filesystem>
#include <vector>
#include <algorithm>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>


crow::json::wvalue VtuData::toJson() const {
    crow::json::wvalue result;

    result["number_of_points"] = numberOfPoints;
    result["number_of_cells"] = numberOfCells;

    if (!points.empty()) {
        std::vector<double> points_vec(points.begin(), points.end());
        result["points"] = std::move(points_vec);
    }

    if (!displacements.empty()) {
        std::vector<double> disp_vec(displacements.begin(), displacements.end());
        result["displacements"] = std::move(disp_vec);
    }

    if (!stresses.empty()) {
        std::vector<double> stress_vec(stresses.begin(), stresses.end());
        result["stresses"] = std::move(stress_vec);
    }

    if (!strains.empty()) {
        std::vector<double> strain_vec(strains.begin(), strains.end());
        result["strains"] = std::move(strain_vec);
    }

    if (!connectivity.empty()) {
        std::vector<int> conn_vec(connectivity.begin(), connectivity.end());
        result["connectivity"] = std::move(conn_vec);
    }

    if (!elementTypes.empty()) {
        std::vector<int> types_vec(elementTypes.begin(), elementTypes.end());
        result["element_types"] = std::move(types_vec);
    }

    if (!displacements.empty() && displacements.size() % 3 == 0) {
        int compCount = 3;
        std::vector<double> dispX, dispY, dispZ;
        dispX.reserve(displacements.size() / 3);
        dispY.reserve(displacements.size() / 3);
        dispZ.reserve(displacements.size() / 3);

        for (size_t i = 0; i + 2 < displacements.size(); i += compCount) {
            dispX.push_back(displacements[i]);
            dispY.push_back(displacements[i + 1]);
            dispZ.push_back(displacements[i + 2]);
        }
        result["disp_x"] = dispX;
        result["disp_y"] = dispY;
        result["disp_z"] = dispZ;
    }

    if (!stresses.empty() && stresses.size() % 6 == 0) {
        int compCount = 6;
        std::vector<double> sx, sy, sz, txy, tyz, tzx;
        sx.reserve(stresses.size() / 6);
        sy.reserve(stresses.size() / 6);

        for (size_t i = 0; i + 5 < stresses.size(); i += compCount) {
            sx.push_back(stresses[i]);
            sy.push_back(stresses[i + 1]);
            sz.push_back(stresses[i + 2]);
            txy.push_back(stresses[i + 3]);
            tyz.push_back(stresses[i + 4]);
            tzx.push_back(stresses[i + 5]);
        }
        result["stress_x"] = sx;
        result["stress_y"] = sy;
        result["stress_z"] = sz;
        result["stress_xy"] = txy;
        result["stress_yz"] = tyz;
        result["stress_zx"] = tzx;
    }

    return result;
}

ApiServer::ApiServer(int port) : port(port) {
    auto& cors = app_.get_middleware<crow::CORSHandler>();
    cors.global()
        .origin("*")
        .headers("Content-Type")
        .methods("GET"_method, "POST"_method, "OPTIONS"_method);
    setupRoutes();
}

VtuData ApiServer::parseVtuFile(const std::string& filename) {
    VtuData data;
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    auto grid = reader->GetOutput();
    if (!grid) {
        std::cerr << "Failed to read file" << std::endl;
        return data;
    }

    auto points = grid->GetPoints();
    if (points) {
        for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            data.points.push_back(p[0]);
            data.points.push_back(p[1]);
            data.points.push_back(p[2]);
        }
        data.numberOfPoints = static_cast<int>(points->GetNumberOfPoints());
    }

    data.numberOfCells = static_cast<int>(grid->GetNumberOfCells());
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
        auto cell = grid->GetCell(i);
        if (!cell) continue;
        auto ids = cell->GetPointIds();
        for (vtkIdType j = 0; j < ids->GetNumberOfIds(); ++j) {
            data.connectivity.push_back(static_cast<int>(ids->GetId(j)));
        }
        data.elementTypes.push_back(static_cast<int>(cell->GetCellType()));
    }

    auto pointData = grid->GetPointData();
    int numArrays = pointData->GetNumberOfArrays();

    vtkDataArray* dispArray = nullptr;
    vtkDataArray* stressArray = nullptr;
    vtkDataArray* strainArray = nullptr;

    for (int i = 0; i < numArrays; ++i) {
        vtkDataArray* arr = pointData->GetArray(i);
        if (!arr) continue;
        std::string name = arr->GetName() ? arr->GetName() : "";

        std::string lname = name;
        std::transform(lname.begin(), lname.end(), lname.begin(), ::tolower);

        if (lname.find("disp") != std::string::npos || lname == "u")
            dispArray = arr;
        else if (lname.find("stress") != std::string::npos || lname == "s")
            stressArray = arr;
        else if (lname.find("strain") != std::string::npos || lname == "e")
            strainArray = arr;
    }

    if (dispArray) {
        int nComp = dispArray->GetNumberOfComponents();
        for (vtkIdType i = 0; i < dispArray->GetNumberOfTuples(); ++i) {
            std::vector<double> tuple(nComp);
            dispArray->GetTuple(i, tuple.data());
            data.displacements.insert(data.displacements.end(), tuple.begin(), tuple.end());
        }
    }

    if (stressArray) {
        int nComp = stressArray->GetNumberOfComponents();
        for (vtkIdType i = 0; i < stressArray->GetNumberOfTuples(); ++i) {
            std::vector<double> tuple(nComp);
            stressArray->GetTuple(i, tuple.data());
            data.stresses.insert(data.stresses.end(), tuple.begin(), tuple.end());
        }
    }

    if (strainArray) {
        int nComp = strainArray->GetNumberOfComponents();
        for (vtkIdType i = 0; i < strainArray->GetNumberOfTuples(); ++i) {
            std::vector<double> tuple(nComp);
            strainArray->GetTuple(i, tuple.data());
            data.strains.insert(data.strains.end(), tuple.begin(), tuple.end());
        }
    }

    return data;
}

std::string ApiServer::getVtuFilename(const std::string& fcFilename) {
    std::string baseName = fcFilename;
    size_t dotPos = baseName.find_last_of('.');
    if (dotPos != std::string::npos && baseName.substr(dotPos) == ".fc") {
        baseName = baseName.substr(0, dotPos);
    }
    return baseName + ".vtu";
}

void ApiServer::setupRoutes() {
    CROW_ROUTE(app_, "/")
        ([]() {
        return "FEM calculate server is running";
        });

    CROW_ROUTE(app_, "/api/analysis").methods("POST"_method)
        ([this](const crow::request& req) {

        try {
            auto json_data = crow::json::load(req.body);
            std::string originalFilename;

            if (json_data && json_data.has("header")) {
                auto header = json_data["header"];
                if (header.has("original_filename")) {
                    originalFilename = header["original_filename"].s();
                }
            }

            this->curr_filename = originalFilename;

            std::ofstream ofs(curr_filename + ".fc");
            ofs << req.body;
            ofs.close();

            logger& log = logger::log();
            log.print("Starting calculation");
            std::cout << "Starting calculation..." << std::endl;

            Calculate solve(originalFilename);
            solve.Solve();

            std::cout << "Calculation completed" << std::endl;
            log.print("Calculation completed");

            std::string resultVtuFile = this->getVtuFilename(originalFilename);

            bool hasResults = !resultVtuFile.empty();

            if (hasResults) {
                VtuData results = this->parseVtuFile(resultVtuFile);
            }

            crow::json::wvalue response;
            response["success"] = true;
            response["message"] = "Analysis completed successfully";
            response["results_available"] = hasResults;

            if (hasResults) {
                response["results_file"] = resultVtuFile;
                response["download_filename"] = this->getVtuFilename(originalFilename);
            }

            return crow::response{ 200, response };
        }
        catch (const std::exception& e) {
            std::cerr << "Error during analysis: " << e.what() << std::endl;
            crow::json::wvalue error;
            error["success"] = false;
            error["message"] = e.what();
            return crow::response{ 400, error };
        }
        });

    CROW_ROUTE(app_, "/api/results")
        ([this](const crow::request& req) {

        try {

            if (curr_filename.empty()) {
                crow::json::wvalue response;
                response["success"] = true;
                response["has_results"] = false;
                response["message"] = "No VTU results file found";
                response["debug_info"] = "No .vtu files found in current directory";
                return crow::response{ 200, response };
            }

            VtuData results = this->parseVtuFile(curr_filename);

            if (results.numberOfPoints == 0) {
                crow::json::wvalue response;
                response["success"] = false;
                response["message"] = "Failed to parse VTU file or file is empty";
                return crow::response{ 500, response };
            }

            crow::json::wvalue response = results.toJson();
            response["success"] = true;
            response["has_results"] = true;
            response["file_path"] = curr_filename;
            response["debug_info"] = "Results parsed successfully";
            response["file_size"] = std::filesystem::file_size(curr_filename);

            return crow::response{ 200, response };

        }
        catch (const std::exception& e) {
            std::cerr << "Error in /api/results: " << e.what() << std::endl;
            crow::json::wvalue error;
            error["success"] = false;
            error["message"] = e.what();
            return crow::response{ 500, error };
        }
        });

    CROW_ROUTE(app_, "/<string>")
        ([](const std::string& filename) {
        if (!std::filesystem::exists(filename)) {
            return crow::response(404, "File not found");
        }
        std::ifstream file(filename, std::ios::binary);
        std::ostringstream oss;
        oss << file.rdbuf();
        crow::response res;
        res.code = 200;
        res.set_header("Content-Type", "application/octet-stream");
        res.write(oss.str());
        res.end();
        return res;
        });

    CROW_ROUTE(app_, "/api/download/vtu")
        ([this](const crow::request& req) {

        try {
            std::string vtuFile = this->getVtuFilename(this->curr_filename);

            if (vtuFile.empty()) {
                crow::json::wvalue error;
                error["success"] = false;
                error["message"] = "No VTU file available for download";
                return crow::response{ 404, error };
            }

            std::ifstream file(vtuFile, std::ios::binary);
            if (!file) {
                crow::json::wvalue error;
                error["success"] = false;
                error["message"] = "Cannot open VTU file";
                return crow::response{ 500, error };
            }

            file.seekg(0, std::ios::end);
            size_t fileSize = file.tellg();
            file.seekg(0, std::ios::beg);

            std::string fileContent(fileSize, '\0');
            file.read(&fileContent[0], fileSize);
            file.close();

            auto response = crow::response{ fileContent };
            response.add_header("Content-Type", "application/xml");
            response.add_header("Content-Disposition",
                "attachment; filename=\"" + vtuFile + "\"");
            response.add_header("Content-Length", std::to_string(fileSize));

            return response;
        }
        catch (const std::exception& e) {
            std::cerr << "Error in /api/download/vtu: " << e.what() << std::endl;
            crow::json::wvalue error;
            error["success"] = false;
            error["message"] = e.what();
            return crow::response{ 500, error };
        }
        });

}


void ApiServer::run() {
    std::cout << "Starting FEM calculate server on port " << port << std::endl;
    std::cout << "Server will be available at: http://localhost:" << port << std::endl;

    try {
        app_.port(port).multithreaded().run();
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: Failed to start server: " << e.what() << std::endl;
        std::cerr << "Make sure port " << port << " is not already in use" << std::endl;
    }
}