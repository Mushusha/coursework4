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
namespace Tensor {
	enum ResType { STRAIN, STRESS, COUNT };
}
namespace Vector {
	enum ResType { DISPLACEMENT, VELOCITY, ACCELERATION, COUNT }; 
}
