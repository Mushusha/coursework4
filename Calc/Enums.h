#pragma once

#include <map>
#include <string>

enum ElemType { TRI = 10, QUAD = 12, TETRA = 1, HEX = 3, WERGE = 6, PYR = 8, INFQUAD = 10};
enum LoadType { PRESSURE = 4, NODEFORCE = 5 };
enum LocVar { KSI, ETA, ZETA };
enum GlobVar { X, Y, Z };

namespace Comp3D {
	enum Comp { XX, YY, ZZ, XY, XZ, YZ };
}
namespace Comp2D {
	enum Comp { XX, YY, XY };
}

enum ResType { STRAIN, STRESS, DISPLACEMENT, VELOCITY, ACCELERATION, COUNT };

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