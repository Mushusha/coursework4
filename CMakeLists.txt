cmake_minimum_required(VERSION 3.15)
set(CMAKE_TOOLCHAIN_FILE "C:/Users/mushusha/vcpkg/scripts/buildsystems/vcpkg.cmake")
project(cw)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "Vcpkg target triplet")
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-3.4.0")
set(VTK_DIR "C:/Users/mushusha/vcpkg/packages/vtk_x64-windows/share/vtk")
set(GLEW_DIR "C:/Users/mushusha/vcpkg/installed/x64-windows/share/glew")

set(VTK_NO_AUTO_FIND_GLEW ON)

find_package(VTK CONFIG REQUIRED)
find_package(GLEW CONFIG REQUIRED)
find_package(pugixml REQUIRED)
find_package(GTest CONFIG REQUIRED)
find_package(Crow CONFIG REQUIRED)

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
Elements/Element.cpp)

source_group(Elements\\Finite FILES
Elements/Finite/Tri.h
Elements/Finite/Tri.cpp
Elements/Finite/Quad.h
Elements/Finite/Quad.cpp
Elements/Finite/Tetra.h
Elements/Finite/Tetra.cpp
Elements/Finite/Hex.h
Elements/Finite/Hex.cpp
Elements/Finite/Wedge.h
Elements/Finite/Wedge.cpp
Elements/Finite/Pyr.h
Elements/Finite/Pyr.cpp)

source_group(Elements\\Infinite FILES
Elements/Infinite/InfQuad.h
Elements/Infinite/InfQuad.cpp
Elements/Infinite/InfHex.h
Elements/Infinite/InfHex.cpp)

source_group(Elements\\Spectral FILES
Elements/Spectral/QuadSEM.h
Elements/Spectral/HexSEM.h
Elements/Spectral/InfQuadSEM.h
Elements/Spectral/InfHexSEM.h)

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

source_group(Server FILES 
Server/Server.h
Server/Server.cpp)

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

source_group(Tests\\Elements FILES
Tests/Elements/TriTests.h
Tests/Elements/TriTests.cpp
Tests/Elements/QuadTests.h
Tests/Elements/QuadTests.cpp
Tests/Elements/HexTests.h
Tests/Elements/HexTests.cpp
Tests/Elements/TetraTests.h
Tests/Elements/TetraTests.cpp
Tests/Elements/WedgeTests.h
Tests/Elements/WedgeTests.cpp
Tests/Elements/PyrTests.h
Tests/Elements/PyrTests.cpp)


source_group(Tests\\Calc FILES
Tests/Calc/MathMVTests.h
Tests/Calc/MathMVTests.cpp
Tests/Calc/DataTests.h
Tests/Calc/DataTests.cpp)

source_group(Tests\\Solve FILES
Tests/Solve/SolveTests.h
Tests/Solve/SolveTests.cpp
Tests/Solve/IntegrationTests.h
Tests/Solve/IntegrationTests.cpp)

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
    
    Elements/Finite/Tri.h
    Elements/Finite/Tri.cpp
    Elements/Finite/Quad.h
    Elements/Finite/Quad.cpp
    Elements/Finite/Tetra.h
    Elements/Finite/Tetra.cpp
    Elements/Finite/Hex.h
    Elements/Finite/Hex.cpp
    Elements/Finite/Wedge.h
    Elements/Finite/Wedge.cpp
    Elements/Finite/Pyr.h
    Elements/Finite/Pyr.cpp
    
    Elements/Infinite/InfQuad.h
    Elements/Infinite/InfQuad.cpp
    Elements/Infinite/InfHex.h
    Elements/Infinite/InfHex.cpp
    
    Elements/Spectral/QuadSEM.h
    Elements/Spectral/HexSEM.h
    Elements/Spectral/InfQuadSEM.h
    Elements/Spectral/InfHexSEM.h

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

    Server/Server.h
    Server/Server.cpp

    Output/Output.h
    Output/Output.cpp
    Output/VTUWriter.h
    Output/VTUWriter.cpp
    Output/PVDWriter.h
    Output/PVDWriter.cpp

    Parser/Parser.h
    Parser/Parser.cpp
	
    Tests/Tests.cpp
    Tests/Tests.h
    Tests/log1.h
    
    Tests/Elements/TriTests.h
    Tests/Elements/TriTests.cpp
    Tests/Elements/QuadTests.h
    Tests/Elements/QuadTests.cpp
    Tests/Elements/HexTests.h
    Tests/Elements/HexTests.cpp
    Tests/Elements/TetraTests.h
    Tests/Elements/TetraTests.cpp
    Tests/Elements/WedgeTests.h
    Tests/Elements/WedgeTests.cpp
    Tests/Elements/PyrTests.h
    Tests/Elements/PyrTests.cpp
    
    Tests/Calc/MathMVTests.h
    Tests/Calc/MathMVTests.cpp
    Tests/Calc/DataTests.h
    Tests/Calc/DataTests.cpp
    
    Tests/Solve/SolveTests.h
    Tests/Solve/SolveTests.cpp
    Tests/Solve/IntegrationTests.cpp
    Tests/Solve/IntegrationTests.h
)

target_link_libraries(Core PRIVATE 
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOXML
    GLEW::GLEW
    pugixml::pugixml
    Crow::Crow
    GTest::gtest
    GTest::gtest_main
)

target_include_directories(Core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Calc
    ${CMAKE_CURRENT_SOURCE_DIR}/Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/Elements/Finite
    ${CMAKE_CURRENT_SOURCE_DIR}/Elements/Infinite
    ${CMAKE_CURRENT_SOURCE_DIR}/Elements/Spectral
    ${CMAKE_CURRENT_SOURCE_DIR}/Node
    ${CMAKE_CURRENT_SOURCE_DIR}/Output
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoothing
    ${CMAKE_CURRENT_SOURCE_DIR}/Server
    ${CMAKE_CURRENT_SOURCE_DIR}/Solve
    ${CMAKE_CURRENT_SOURCE_DIR}/Parser
    ${CMAKE_CURRENT_SOURCE_DIR}/Tests
    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Node
    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Calc
    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Solve
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)
