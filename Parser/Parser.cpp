#include "Parser.h"

inline bool is_base64(unsigned char c) {
	return (isalnum(c) || (c == '+') || (c == '/'));
}
void base64_decode(std::string& encoded_string)
{
	std::string decoded_string;
	int in_len = static_cast<int>(encoded_string.size());
	int i = 0;
	int j = 0;
	int in_ = 0;
	unsigned char char_array_4[4], char_array_3[3];
	decoded_string.clear();

	while (in_len-- && (encoded_string[in_] != '=') && is_base64(encoded_string[in_])) {
		char_array_4[i++] = encoded_string[in_]; in_++;
		if (i == 4) {
			for (i = 0; i < 4; i++)
				char_array_4[i] = static_cast<unsigned char>(base64_chars.find(char_array_4[i]));

			char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
			char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
			char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

			for (i = 0; (i < 3); i++)
				decoded_string += char_array_3[i];
			i = 0;
		}
	}

	if (i) {
		for (j = i; j < 4; j++)
			char_array_4[j] = 0;

		for (j = 0; j < 4; j++)
			char_array_4[j] = static_cast<unsigned char>(base64_chars.find(char_array_4[j]));

		char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
		char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
		char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

		for (j = 0; (j < i - 1); j++) decoded_string += char_array_3[j];
	}
	swap(encoded_string, decoded_string);
}

lint ReadInt(string& data, size_t i) {
	return *reinterpret_cast<const int*>(data.c_str() + (i) * sizeof(int));
}
int ReadUInt8_t(string& data, size_t i) {
	return *reinterpret_cast<const int*>(data.c_str() + (i) * sizeof(uint8_t));
}
double ReadDouble(string& data, size_t i) {
	return *reinterpret_cast<const double*>(data.c_str() + (i) * sizeof(double));
}
bool ReadBool(string& data, size_t i) {
	return *reinterpret_cast<const bool*>(data.c_str() + (i) * sizeof(bool));
}


void block::read(json block) {
	this->id = block["id"];
	this->cs_id = block["cs_id"];
	this->mat_id = block["material_id"];
	this->prop_id = block["property_id"];
}
void coordinate::read(json cs) {
	this->id = cs["id"];

	if (cs["type"] == "cartesian") {
		this->type = type_cs::Cartesian;
	}
	else if (cs["type"] == "cylindrical") {
		this->type = type_cs::Cylindrical;
	}
	else {
		this->type = type_cs::Unknown;
		std::cout << "ERROR: unknown coordinate type";
	}

	std::string dir1_s = cs["dir1"];
	std::string dir2_s = cs["dir2"];
	std::string origin_s = cs["origin"];
	base64_decode(dir1_s);
	base64_decode(dir2_s);
	base64_decode(origin_s);
	for (uint8_t i = 0; i < 3; i++) {
		this->dir1[i] = ReadDouble(dir1_s, i);
		this->dir2[i] = ReadDouble(dir2_s, i);
		this->origin[i] = ReadDouble(origin_s, i);
	}
}
void load::read(json load) {
	this->id = load["id"];
	//this->cs = load["cs"];
	this->type = load["type"];
	this->size = load["apply_to_size"];
	if (load.contains("inf"))
		this->inf = load["inf"];

	this->apply_to.resize(2 * this->size);
	std::string apply_to_s = load["apply_to"];
	base64_decode(apply_to_s);
	for (size_t i = 0; i < this->apply_to.size(); i++)
		this->apply_to[i] = ReadInt(apply_to_s, i);

	for (uint8_t i = 0; i < load["data"].size(); i++) {
		std::string s_data = load["data"][i];
		base64_decode(s_data);
		this->data[i] = ReadDouble(s_data, 0);
	}
}

void nodesets::read(json nodesets) {
	this->id = nodesets["id"];
	this->size = nodesets["apply_to_size"];
	this->apply_to.resize(this->size);
	std::string apply_to_s = nodesets["apply_to"];
	base64_decode(apply_to_s);
	for (size_t i = 0; i < this->apply_to.size(); i++)
		this->apply_to[i] = ReadInt(apply_to_s, i);
}

void sidesets::read(json sidesets) {
	this->id = sidesets["id"];
	this->size = sidesets["apply_to_size"];
	this->apply_to.resize(2 * this->size);
	this->load = sidesets["load"];
	std::string apply_to_s = sidesets["apply_to"];
	base64_decode(apply_to_s);
	for (size_t i = 0; i < this->apply_to.size(); i++)
		this->apply_to[i] = ReadInt(apply_to_s, i);
}

void material::read(json mat) {
	this->id = mat["id"];
	this->type = mat["elasticity"][0]["type"];

	for (uint8_t i = 0; i < mat["elasticity"][0]["constants"].size(); i++) { // 1 - E, 2 - nu, 3 - rho
		std::string tmp = mat["elasticity"][0]["constants"][i];
		base64_decode(tmp);
		this->constants[i] = ReadDouble(tmp, 0);
	}
	if (mat.contains("common")) {
		std::string tmp = mat["common"][0]["constants"][0];
		base64_decode(tmp);
		this->constants[2] = ReadDouble(tmp, 0);
	}
}
void mesh::read(json mesh) {
	this->elems_count = mesh["elems_count"];

	{
		std::string s_elem_blocks = mesh["elem_blocks"];
		std::string s_elem_orders = mesh["elem_orders"];
		std::string s_elem_parent_ids = mesh["elem_parent_ids"];
		std::string s_elem_types = mesh["elem_types"];
		std::string s_elemids = mesh["elemids"];
		std::string s_elems = mesh["elems"];
		base64_decode(s_elem_blocks);
		base64_decode(s_elem_orders);
		base64_decode(s_elem_parent_ids);
		base64_decode(s_elem_types);
		base64_decode(s_elemids);
		base64_decode(s_elems);
		this->elem_blocks.resize(this->elems_count);
		this->elem_orders.resize(this->elems_count);
		this->elem_parent_ids.resize(this->elems_count);
		this->elem_types.resize(this->elems_count);
		this->elem_id.resize(this->elems_count);
		lint count = 0;

		for (lint i = 0; i < this->elems_count; i++) {
			this->elem_blocks[i] = static_cast<int>(ReadInt(s_elem_blocks, i));
			this->elem_orders[i] = static_cast<uint8_t>(ReadUInt8_t(s_elem_orders, i));
			this->elem_parent_ids[i] = static_cast<int>(ReadInt(s_elem_parent_ids, i));
			this->elem_types[i] = static_cast<uint8_t>(ReadUInt8_t(s_elem_types, i));
			this->elem_id[i] = ReadInt(s_elemids, i);
			count += count_nodes(s_elem_types[i]);
		}

		this->elem_nodes.resize(count);
		for (lint i = 0; i < count; i++) {
			this->elem_nodes[i] = ReadInt(s_elems, i);
		}

		this->nodes_count = mesh["nodes_count"];
		std::string s_nids = mesh["nids"];
		std::string s_nodes = mesh["nodes"];
		base64_decode(s_nids);
		base64_decode(s_nodes);
		this->node_id.resize(this->nodes_count);
		this->nodes_coord.resize(this->nodes_count * 3);
		for (lint i = 0; i < this->nodes_count; i++) {
			this->node_id[i] = static_cast<lint>(ReadInt(s_nids, i));
			this->nodes_coord[i * 3] = (ReadDouble(s_nodes, i * 3));
			this->nodes_coord[i * 3 + 1] = (ReadDouble(s_nodes, i * 3 + 1));
			this->nodes_coord[i * 3 + 2] = (ReadDouble(s_nodes, i * 3 + 2));
		}
	}
}
uint8_t mesh::count_nodes(uint8_t elem_t) {
	switch (elem_t)
	{
	case 10:
		return 3;
		break;
	default:
		return 4;
	}
}

void restraints::read(json restraints) {
	this->id = restraints["id"];
	this->cs = restraints["cs"];
	this->size = restraints["apply_to_size"];

	this->apply_to.resize(this->size);
	std::string s_apply_to = restraints["apply_to"];
	base64_decode(s_apply_to);
	for (size_t i = 0; i < this->apply_to.size(); i++)
		this->apply_to[i] = ReadInt(s_apply_to, i);

	for (size_t i = 0; i < restraints["data"].size(); i++) {
		std::string s_data = restraints["data"][i];
		base64_decode(s_data);
		this->data[i] = ReadDouble(s_data, i);
	}

	for (size_t i = 0; i < restraints["flag"].size(); i++)
		this->flag[i] = restraints["flag"][i];
}
void settings::read(json block) {
	this->dimensions = block["dimensions"];
	if (this->dimensions == "2D" && block.contains("plane_state"))
		this->plane_state = block["plane_state"];
	this->analisys_type = block["type"];
	if (this->analisys_type == "dynamic") {
		this->max_time = block["dynamics"]["max_time"];
		this->d = block["damping"]["mass_matrix"];
	}
}

void Parser::read(string name) {
	this->filename = name;
	std::ifstream fc_file(this->filename, std::ios::in);
	if (!fc_file)
		throw std::runtime_error("cannot open fc file: " + this->filename);
	auto _root = nlohmann::json::parse(fc_file);
	fc_file.close();

	this->block.resize(_root["blocks"].size());
	for (int i = 0; i < this->block.size(); i++) {
		this->block[i].read(_root["blocks"][i]);
	}
	
	this->coordinate.resize(_root["coordinate_systems"].size());
	for (auto i = 0; i < this->coordinate.size(); i++) {
		this->coordinate[i].read(_root["coordinate_systems"][i]);
	}
	
	this->load.resize(_root["loads"].size());
	for (int i = 0; i < this->load.size(); i++) {
		this->load[i].read(_root["loads"][i]);
	}
	
	this->nodesets.resize(_root["sets"]["nodesets"].size());
	for (int i = 0; i < this->nodesets.size(); i++) {
		this->nodesets[i].read(_root["sets"]["nodesets"][i]);
	}
	
	this->sidesets.resize(_root["sets"]["sidesets"].size());
	for (int i = 0; i < this->sidesets.size(); i++) {
		this->sidesets[i].read(_root["sets"]["sidesets"][i]);
	}
	
	this->material.resize(_root["materials"].size());
	for (int i = 0; i < this->material.size(); i++) {
		this->material[i].read(_root["materials"][i]);
	}
	
	this->mesh.read(_root["mesh"]);
	
	this->restraints.resize(_root["restraints"].size());
	for (int i = 0; i < this->restraints.size(); i++) {
		this->restraints[i].read(_root["restraints"][i]);
	}
	
	this->settings.read(_root["settings"]);
}