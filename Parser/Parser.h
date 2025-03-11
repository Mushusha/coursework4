#pragma once
#include "json.hpp"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using json = nlohmann::json;
typedef long long int  lint;

//Decode
inline static const std::string base64_chars =
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz"
"0123456789+/";
inline bool is_base64(unsigned char c);
void base64_decode(std::string& encoded_string);

//Read type data
inline double ReadDouble(string& data, size_t i);
inline lint ReadInt(string& data, size_t i);
inline int ReadUInt8_t(string& data, size_t i);
inline bool ReadBool(string& data, size_t i);

enum type_cs {
	Cartesian = 0,
	Cylindrical = 1,
	Unknown = 2
};

struct block {
	int id;
	int cs_id;
	int mat_id;
	int prop_id;
	void read(json block);
};

struct coordinate {
	int id;
	double dir1[3];
	double dir2[3];
	double origin[3];
	type_cs type;
	void read(json cs);
};

struct load {
	int id;
	std::vector<lint> apply_to;
	size_t size;
	int cs;
	std::array<double, 6> data;
	int type;
	bool inf = false;
	void read(json load);
};

struct nodesets {
	int id;
	std::vector<lint> apply_to;
	size_t size;
	void read(json nodesets);
};

struct sidesets {
	int id;
	std::vector<lint> apply_to;
	size_t size;
	int load;
	void read(json sidesets);
};

struct material {
	int id;
	int type;
	std::array<double, 3> constants;
	void read(json material);
};

struct mesh {
	lint elems_count;
	std::vector<int> elem_blocks;
	std::vector<uint8_t> elem_orders;
	std::vector<int> elem_parent_ids;
	std::vector<uint8_t> elem_types;
	std::vector<lint> elem_id;
	std::vector<lint> elem_nodes;
	lint nodes_count;
	std::vector<lint> node_id;
	std::vector<double> nodes_coord;

	inline uint8_t count_nodes(uint8_t elem_types);
	void read(json mesh);
};

struct restraints {
	int id;
	std::vector<lint> apply_to;
	int size;
	int cs;
	double data[6];
	int flag[6];
	void read(json restraints);
};

struct settings {
	std::string plane_state;
	std::string dimensions;
	void read(json settings);
};

class Parser {
public:
	Parser() = default;
	std::vector<block> block;
	std::vector<coordinate> coordinate;
	std::vector<load> load;
	std::vector<material> material;
	std::vector<nodesets> nodesets;
	std::vector<sidesets> sidesets;
	mesh mesh;
	std::vector<restraints> restraints;
	settings settings;
	void read(string filename);

private:
	string filename;
};
#pragma once
