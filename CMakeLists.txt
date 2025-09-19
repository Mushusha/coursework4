cmake_minimum_required(VERSION 3.12)
set(CMAKE_TOOLCHAIN_FILE "C:/Users/mushusha/vcpkg/scripts/buildsystems/vcpkg.cmake")
project(cw)

set(GLEW_DIR "C:/Users/mushusha/vcpkg/installed/x64-windows/share/glew")

set(VTK_NO_AUTO_FIND_GLEW ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "Vcpkg target triplet")
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-3.4.0")
set(VTK_DIR "C:/Users/mushusha/vcpkg/packages/vtk_x64-windows/share/vtk")

find_package(VTK CONFIG REQUIRED)
find_package(GLEW CONFIG REQUIRED)
find_package(pugixml REQUIRED)
find_package(GTest REQUIRED) 

source_group(Calc FILES 
Calc/main.cpp
Calc/Calc.h
Calc/Calc.cpp
Calc/Data.h
Calc/Data.cpp
Calc/MathMV.h
Calc/MathMV.cpp
Calc/Enums.h)

source_group(Node FILES
Node/Node.h
Node/Node.cpp)

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
Smoothing/Interpolation.cpp)

source_group(Output FILES 
Output/Output.h
Output/Output.cpp
Output/VTUWriter.h
Output/VTUWriter.cpp
Output/PVDWriter.h
Output/PVDWriter.cpp)

source_group(Parser FILES 
Parser/Parser.h
Parser/Parser.cpp)

source_group(Tests FILES 
Tests/Tests.h
Tests/Tests.cpp
Tests/log1.h)

add_executable(Core
    Calc/main.cpp
    Calc/Calc.h
    Calc/Calc.cpp
    Calc/Data.h
    Calc/Data.cpp
    Calc/MathMV.h
    Calc/MathMV.cpp
    Calc/Enums.h

    Node/Node.h
    Node/Node.cpp

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
    Elements/Hex.cpp

    Solve/Solver.h
    Solve/Solver.cpp
    Solve/Statics.h
    Solve/Statics.cpp
    Solve/Dynamics.h
    Solve/Dynamics.cpp

    Smoothing/Smoothing.h
    Smoothing/Smoothing.cpp
    Smoothing/Interpolation.h
    Smoothing/Interpolation.cpp

    Output/Output.h
    Output/Output.cpp
    Output/VTUWriter.h
    Output/VTUWriter.cpp
    Output/PVDWriter.h
    Output/PVDWriter.cpp

    Parser/Parser.h
    Parser/Parser.cpp

    Tests/Tests.h
    Tests/Tests.cpp
    Tests/log1.h
)

target_link_libraries(Core PRIVATE 
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOXML
    GLEW::GLEW
    pugixml::pugixml
    GTest::GTest
    GTest::Main
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