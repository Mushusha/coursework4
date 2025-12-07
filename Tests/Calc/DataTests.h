#pragma once

#include <gtest/gtest.h>
#include <memory>

#include "Data.h"
#include "Parser.h"
#include "Node/Node.h"
#include "Tri.h"
#include "Quad.h"
#include "Hex.h"
#include "Tetra.h"
#include "Wedge.h"
#include "Pyr.h"

class MockParser : public Parser {
public:
    MockParser() {
        settings.dimensions = "2D";
        settings.plane_state = "p-stress";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 4;
        mesh.node_id = { 1, 2, 3, 4 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { QUAD };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParser3D : public Parser {
public:
    MockParser3D() {
        settings.dimensions = "3D";
        settings.plane_state = "";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 8;
        mesh.node_id = { 1, 2, 3, 4, 5, 6, 7, 8 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            0.0, 1.0, 1.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { HEX };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4, 5, 6, 7, 8 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParserTri : public Parser {
public:
    MockParserTri() {
        settings.dimensions = "2D";
        settings.plane_state = "p-stress";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 3;
        mesh.node_id = { 1, 2, 3 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { TRI };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParserTetra : public Parser {
public:
    MockParserTetra() {
        settings.dimensions = "3D";
        settings.plane_state = "";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 4;
        mesh.node_id = { 1, 2, 3, 4 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { TETRA };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParserQuad : public Parser {
public:
    MockParserQuad() {
        settings.dimensions = "2D";
        settings.plane_state = "p-stress";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 4;
        mesh.node_id = { 1, 2, 3, 4 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { QUAD };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParserHex : public Parser {
public:
    MockParserHex() {
        settings.dimensions = "3D";
        settings.plane_state = "";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 8;
        mesh.node_id = { 1, 2, 3, 4, 5, 6, 7, 8 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            0.0, 1.0, 1.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { HEX };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4, 5, 6, 7, 8 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParserWedge : public Parser {
public:
    MockParserWedge() {
        settings.dimensions = "3D";
        settings.plane_state = "";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 6;
        mesh.node_id = { 1, 2, 3, 4, 5, 6 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            0.0, 1.0, 1.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { WEDGE };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4, 5, 6 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParserPyr : public Parser {
public:
    MockParserPyr() {
        settings.dimensions = "3D";
        settings.plane_state = "";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 5;
        mesh.node_id = { 1, 2, 3, 4, 5 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.5, 0.5, 1.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { PYR };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4, 5 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        this->restraints.clear();
        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class MockParserWithConstraints : public Parser {
public:
    MockParserWithConstraints() {
        settings.dimensions = "2D";
        settings.plane_state = "p-stress";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 4;
        mesh.node_id = { 1, 2, 3, 4 };
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { QUAD };
        mesh.elem_orders = { 1 };
        mesh.elem_nodes = { 1, 2, 3, 4 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        ::restraints rest;
        rest.id = 1;
        rest.apply_to = { 1, 4 };
        rest.size = 2;
        rest.cs = 0;
        rest.data[0] = 0.0;
        rest.data[1] = 0.0;
        rest.data[2] = 0.0;
        rest.flag[0] = 1;
        rest.flag[1] = 1;
        rest.flag[2] = 0;
        rest.flag[3] = 0;
        rest.flag[4] = 0;
        rest.flag[5] = 0;
        this->restraints.push_back(rest);

        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

class DataTest : public ::testing::Test {
protected:
    void SetUp() override {
        parser2D = std::make_shared<MockParser>();
        parser3D = std::make_shared<MockParser3D>();
    }

    std::shared_ptr<MockParser> parser2D;
    std::shared_ptr<MockParser3D> parser3D;
};

class DataCopyMoveTest : public ::testing::Test {
protected:
    void SetUp() override {
        parser = std::make_shared<MockParser>();
    }

    std::shared_ptr<MockParser> parser;
};

class DataConstraintsTest : public ::testing::Test {
protected:
    void SetUp() override {
        parser = std::make_shared<MockParserWithConstraints>();
    }

    std::shared_ptr<MockParserWithConstraints> parser;
};

class DataElementTypesTest : public ::testing::Test {
protected:
    void SetUp() override {
        parserTri = std::make_shared<MockParserTri>();
        parserQuad = std::make_shared<MockParserQuad>();
        parserTetra = std::make_shared<MockParserTetra>();
        parserHex = std::make_shared<MockParserHex>();
        parserWedge = std::make_shared<MockParserWedge>();
        parserPyr = std::make_shared<MockParserPyr>();
    }

    std::shared_ptr<MockParserTri> parserTri;
    std::shared_ptr<MockParserQuad> parserQuad;
    std::shared_ptr<MockParserTetra> parserTetra;
    std::shared_ptr<MockParserHex> parserHex;
    std::shared_ptr<MockParserWedge> parserWedge;
    std::shared_ptr<MockParserPyr> parserPyr;
};
