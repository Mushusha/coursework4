cmake_minimum_required(VERSION 3.10)
project(cw)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

source_group(Calc FILES 
Calc/main.cpp
Calc/Data.h
Calc/Data.cpp
Calc/MathMV.h
Calc/MathMV.cpp
Calc/Node.h
Calc/Node.cpp
Calc/Enums.h)

source_group(Elements FILES 
Elements/Element.h
Elements/Element.cpp
Elements/Tri.h
Elements/Tri.cpp
Elements/Quad.h
Elements/Quad.cpp
Elements/InfQuad.cpp
Elements/Tetra.h
Elements/Tetra.cpp
Elements/Hex.h
Elements/Hex.cpp)

source_group(Solve FILES 
Solve/Solver.h
Solve/Solver.cpp
Solve/Statics.h
Solve/Statics.cpp
Solve/Dynamics.h
Solve/Dynamics.cpp
Solve/Smoothing.h
Solve/Smoothing.cpp
Solve/Fabric.h)

source_group(Smoothing FILES 
Smoothing/Smoothing.h
Smoothing/Smoothing.cpp
Smoothing/Interpolation.h
Smoothing/Interpolation.cpp
Smoothing/Fabric.h)

source_group(Parser FILES 
Parser/Parser.h
Parser/Parser.cpp)


source_group(Tests FILES 
Tests/Tests.h
Tests/Tests.cpp
Tests/log1.h)


set(EIGEN3_INCLUDE_DIR "eigen-3.4.0")

add_executable(Core
    Calc/main.cpp
    Calc/Data.h
    Calc/Data.cpp
    Calc/MathMV.h
    Calc/MathMV.cpp
    Calc/Node.h
    Calc/Node.cpp
    Calc/Enums.h

    Elements/Element.h
    Elements/Element.cpp
    Elements/Tri.h
    Elements/Tri.cpp
    Elements/Quad.cpp
    Elements/InfQuad.cpp
    Elements/Tetra.h
    Elements/Tetra.cpp
    Elements/Hex.h
    Elements/Hex.cpp

    Solve/Solver.h
    Solve/Solver.cpp
    Solve/Statics.h
    Solve/Statics.cpp
    Solve/Dynamics.h
    Solve/Dynamics.cpp
    Solve/Fabric.h

    Smoothing/Smoothing.h
    Smoothing/Smoothing.cpp
    Smoothing/Interpolation.h
    Smoothing/Interpolation.cpp
    Smoothing/Fabric.h

    Parser/Parser.h
    Parser/Parser.cpp

    Tests/Tests.h
    Tests/Tests.cpp
    Tests/log1.h
)

target_include_directories(Core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Calc
    ${CMAKE_CURRENT_SOURCE_DIR}/Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/Solve
    ${CMAKE_CURRENT_SOURCE_DIR}/Parser
    ${CMAKE_CURRENT_SOURCE_DIR}/Tests
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)

