#include "PVDWriter.h"

std::vector<TimestepFile> parse_timesteps(const std::vector<fs::path>& files) {
    std::vector<TimestepFile> result;

    for (const auto& file : files) {
        std::string stem = file.stem().string();
        size_t last_underscore = stem.find_last_of('_');

        if (last_underscore == std::string::npos) continue;

        std::string time_str = stem.substr(last_underscore + 1);

        try {
            double time = std::stod(time_str);
            result.push_back({ file.filename().string(), time });
        }
        catch (...) {
            continue;
        }
    }

    return result;
}

void create_pvd_file(const std::vector<TimestepFile>& timesteps, const fs::path& output_path) {
    pugi::xml_document doc;

    pugi::xml_node vtk_file = doc.append_child("VTKFile");
    vtk_file.append_attribute("type") = "Collection";
    vtk_file.append_attribute("version") = "0.1";

    pugi::xml_node collection = vtk_file.append_child("Collection");

    for (const auto& ts : timesteps) {
        pugi::xml_node dataset = collection.append_child("DataSet");
        dataset.append_attribute("timestep") = ts.time;
        dataset.append_attribute("part") = "0";
        dataset.append_attribute("file") = ts.filename.c_str();
    }

    if (!doc.save_file(output_path.string().c_str())) {
        throw std::runtime_error("Failed to save PVD file");
    }
}
