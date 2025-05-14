#pragma once
#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <pugixml.hpp>

namespace fs = std::filesystem;

struct TimestepFile {
    std::string filename;
    double time;
};

std::vector<TimestepFile> parse_timesteps(const std::vector<fs::path>& files);

void create_pvd_file(const std::vector<TimestepFile>& timesteps, const fs::path& output_path);