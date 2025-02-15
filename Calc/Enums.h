#pragma once

#include <map>
#include <string>

enum ElemType { TRI = 10, QUAD = 12, TETRA = 1, HEX = 3, WERGE = 6, PYR = 8, INFQUAD = 10};
enum LoadType { PRESSURE = 4 };
enum LocVar { KSI, ETA, ZETA };
enum GlobVar { X, Y, Z };
enum Comp3D { XX, YY, ZZ, XY, XZ, YZ };
enum Comp2D { XX_2D, YY_2D, XY_2D };
enum ResType { DISPLACEMENT, STRAIN, STRESS, COUNT };
