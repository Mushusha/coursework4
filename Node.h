#pragma once

class Node {
private:
	int id;
	double x;
	double y;
	double z;
public:
	Node() {};
	Node(int id, std::array <double, 3> coords) {
		this->id = id;
		this->x = coords[0];
		this->y = coords[1];
		this->z = coords[2];
	};
	 virtual ~Node() = default;

	std::map <int, double> constrains;

	int getID() { return id; }
	double getX() { return x; }
	double getY() { return y; }
	double getZ() { return z; }
};