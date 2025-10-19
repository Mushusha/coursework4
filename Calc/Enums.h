#pragma once

#include <map>
#include <string>

enum ElemType { TRI = 10, QUAD = 12, TETRA = 1, HEX = 3, WEDGE = 6, PYR = 8, QUADSEM = 13, INFQUAD = 14 };
enum LoadType { PRESSURE = 3, NODEFORCE = 5, BERLAGE = 20 };
enum LocVar { KSI, ETA, ZETA };
enum GlobVar { X, Y, Z };

namespace Comp3D {
	enum Comp { XX, YY, ZZ, XY, XZ, YZ };
}
namespace Comp2D {
	enum Comp { XX, YY, XY };
}

enum ResType { STRAIN, STRESS, DISPLACEMENT, VELOCITY, ACCELERATION, COUNT };

inline constexpr int count_nodes(ElemType elem_t) {
	switch (elem_t) {
	case TRI:
		return 3;
		break;
	case QUAD:
		return 4;
		break;
	case TETRA:
		return 4;
		break;
	case HEX:
		return 8;
		break;
	case WEDGE:
		return 6;
		break;
	case PYR:
		return 5;
		break;
	case QUADSEM:
		return 8;
		break;
	default:
		throw std::runtime_error("Error: wrong element type");
	}
}

inline bool isTensor(ResType type) {
	if (type == STRAIN ||
		type == STRESS)
		return true;
	return false;
}
inline bool isVector(ResType type) {
	if (type == DISPLACEMENT ||
		type == VELOCITY ||
		type == ACCELERATION)
		return true;
	return false;
}
inline int numComp(ResType type, int dim) {
	switch (type) {
		case STRAIN:
		case STRESS:
			if (dim == 2)
				return 3;
			if (dim == 3)
				return 6;
		case DISPLACEMENT:
		case VELOCITY:
		case ACCELERATION:
			if (dim == 2)
				return 2;
			if (dim == 3)
				return 3;
		case COUNT:
		default:
			return 0;
	}
}
inline constexpr int numTensor() {
	return 2;
}
inline constexpr int numVector() {
	return 3;
}