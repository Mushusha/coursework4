// #define DISABLE_DATA_TESTS

#ifndef DISABLE_DATA_TESTS

#include "DataTests.h"
#include "InfQuad.h"
#include "InfHex.h"

TEST_F(DataTest, Constructor2D_CreatesCorrectDimension) {

    Data data(parser2D);
    
    EXPECT_EQ(data.dim, 2);
}

TEST_F(DataTest, Constructor3D_CreatesCorrectDimension) {

    Data data(parser3D);
    
    EXPECT_EQ(data.dim, 3);
}

TEST_F(DataTest, Constructor2D_CreatesCorrectNodeCount) {

    Data data(parser2D);
    
    EXPECT_EQ(data.nodes_count(), 4);
}

TEST_F(DataTest, Constructor3D_CreatesCorrectNodeCount) {

    Data data(parser3D);
    
    EXPECT_EQ(data.nodes_count(), 8);
}

TEST_F(DataTest, Constructor2D_CreatesCorrectElementCount) {

    Data data(parser2D);
    
    EXPECT_EQ(data.elements_count(), 1);
}

TEST_F(DataTest, Constructor3D_CreatesCorrectElementCount) {

    Data data(parser3D);
    
    EXPECT_EQ(data.elements_count(), 1);
}

TEST_F(DataTest, Constructor2D_SetsAnalysisType) {

    Data data(parser2D);
    
    EXPECT_EQ(data.analisys_type, "static");
}

TEST_F(DataTest, Constructor2D_SetsSettingsFromParser) {

    Data data(parser2D);
    
    EXPECT_DOUBLE_EQ(data.damping, 0.0);
    EXPECT_DOUBLE_EQ(data.max_time, 1.0);
    EXPECT_EQ(data.max_iter, 100);
    EXPECT_EQ(data.iter_res_output, 1);
}

TEST_F(DataTest, GetNode_ReturnsCorrectNode) {

    Data data(parser2D);
    
    auto node = data.get_node(0);
    ASSERT_NE(node, nullptr);
    EXPECT_EQ(node->getID(), 1);
    EXPECT_DOUBLE_EQ(node->getX(), 0.0);
    EXPECT_DOUBLE_EQ(node->getY(), 0.0);
}

TEST_F(DataTest, GetNode_ReturnsCorrectCoordinates) {

    Data data(parser2D);
    
    auto node2 = data.get_node(1);
    EXPECT_DOUBLE_EQ(node2->getX(), 1.0);
    EXPECT_DOUBLE_EQ(node2->getY(), 0.0);
    
    auto node3 = data.get_node(2);
    EXPECT_DOUBLE_EQ(node3->getX(), 1.0);
    EXPECT_DOUBLE_EQ(node3->getY(), 1.0);
}

TEST_F(DataTest, GetElem_ReturnsCorrectElement) {

    Data data(parser2D);
    
    auto elem = data.get_elem(0);
    ASSERT_NE(elem, nullptr);
    EXPECT_EQ(elem->get_type(), ElemType::QUAD);
}

TEST_F(DataTest, GetElem3D_ReturnsCorrectElement) {

    Data data(parser3D);
    
    auto elem = data.get_elem(0);
    ASSERT_NE(elem, nullptr);
    EXPECT_EQ(elem->get_type(), ElemType::HEX);
}

TEST_F(DataTest, GetElements_ReturnsAllElements) {

    Data data(parser2D);
    
    auto elements = data.get_elements();
    EXPECT_EQ(elements.size(), 1);
}

TEST_F(DataTest, GetNodes_ReturnsAllNodes) {

    Data data(parser2D);
    
    auto nodes = data.get_nodes();
    EXPECT_EQ(nodes.size(), 4);
}

TEST_F(DataCopyMoveTest, CopyConstructor_CopiesDim) {

    Data original(parser);
    Data copy(original);
    
    EXPECT_EQ(copy.dim, original.dim);
}

TEST_F(DataCopyMoveTest, CopyConstructor_CopiesNodeCount) {

    Data original(parser);
    Data copy(original);
    
    EXPECT_EQ(copy.nodes_count(), original.nodes_count());
}

TEST_F(DataCopyMoveTest, CopyConstructor_CopiesElementCount) {

    Data original(parser);
    Data copy(original);
    
    EXPECT_EQ(copy.elements_count(), original.elements_count());
}

TEST_F(DataCopyMoveTest, CopyConstructor_CopiesSettings) {

    Data original(parser);
    Data copy(original);
    
    EXPECT_EQ(copy.analisys_type, original.analisys_type);
    EXPECT_DOUBLE_EQ(copy.damping, original.damping);
    EXPECT_DOUBLE_EQ(copy.max_time, original.max_time);
    EXPECT_EQ(copy.max_iter, original.max_iter);
}

TEST_F(DataCopyMoveTest, CopyConstructor_IndependentCopy) {

    Data original(parser);
    Data copy(original);
    
    copy.dim = 3;
    EXPECT_EQ(original.dim, 2);
}

TEST_F(DataCopyMoveTest, CopyAssignment_CopiesCorrectly) {

    Data original(parser);
    Data copy(parser);
    
    copy = original;
    
    EXPECT_EQ(copy.dim, original.dim);
    EXPECT_EQ(copy.nodes_count(), original.nodes_count());
    EXPECT_EQ(copy.elements_count(), original.elements_count());
}

TEST_F(DataCopyMoveTest, CopyAssignment_SelfAssignment) {

    Data data(parser);
    int original_dim = data.dim;
    
    data = data;
    
    EXPECT_EQ(data.dim, original_dim);
}

TEST_F(DataCopyMoveTest, MoveConstructor_MovesData) {

    Data original(parser);
    int original_dim = original.dim;
    int original_node_count = original.nodes_count();
    
    Data moved(std::move(original));
    
    EXPECT_EQ(moved.dim, original_dim);
    EXPECT_EQ(moved.nodes_count(), original_node_count);
}

TEST_F(DataCopyMoveTest, MoveConstructor_SourceIsCleared) {

    Data original(parser);
    
    Data moved(std::move(original));
    
    EXPECT_EQ(original.dim, 0);
    EXPECT_EQ(original.nodes_count(), 0);
    EXPECT_EQ(original.elements_count(), 0);
}

TEST_F(DataCopyMoveTest, MoveAssignment_MovesData) {

    Data original(parser);
    Data target(parser);
    
    int original_dim = original.dim;
    
    target = std::move(original);
    
    EXPECT_EQ(target.dim, original_dim);
}

TEST_F(DataCopyMoveTest, MoveAssignment_SelfAssignment) {

    Data data(parser);
    int original_dim = data.dim;
    
    data = std::move(data);
    
    EXPECT_EQ(data.dim, original_dim);
}

TEST_F(DataConstraintsTest, ConstraintsAppliedToNodes) {

    Data data(parser);
    
    auto node1 = data.get_node(0);
    EXPECT_FALSE(node1->constraints.empty());
    EXPECT_EQ(node1->constraints.size(), 2);
}

TEST_F(DataConstraintsTest, ConstraintsHaveCorrectValues) {

    Data data(parser);
    
    auto node1 = data.get_node(0);
    EXPECT_DOUBLE_EQ(node1->constraints[0], 0.0);
    EXPECT_DOUBLE_EQ(node1->constraints[1], 0.0);
}

TEST_F(DataConstraintsTest, UnconstrainedNodesAreEmpty) {

    Data data(parser);
    
    auto node2 = data.get_node(1);
    EXPECT_TRUE(node2->constraints.empty());
}

TEST_F(DataElementTypesTest, TriElement_CreatedCorrectly) {

    Data data(parserTri);
    
    EXPECT_EQ(data.elements_count(), 1);
    EXPECT_EQ(data.get_elem(0)->get_type(), ElemType::TRI);
    EXPECT_EQ(data.get_elem(0)->nodes_count(), 3);
}

TEST_F(DataElementTypesTest, TetraElement_CreatedCorrectly) {

    Data data(parserTetra);
    
    EXPECT_EQ(data.elements_count(), 1);
    EXPECT_EQ(data.get_elem(0)->get_type(), ElemType::TETRA);
    EXPECT_EQ(data.get_elem(0)->nodes_count(), 4);
}

TEST_F(DataTest, ElementHasCorrectMaterialConstants) {

    Data data(parser2D);
    
    auto elem = data.get_elem(0);
    EXPECT_DOUBLE_EQ(elem->get_E(), 210e9);
    EXPECT_DOUBLE_EQ(elem->get_nu(), 0.3);
    EXPECT_DOUBLE_EQ(elem->get_rho(), 7850.0);
}

TEST_F(DataTest, Element3DHasCorrectMaterialConstants) {

    Data data(parser3D);
    
    auto elem = data.get_elem(0);
    EXPECT_DOUBLE_EQ(elem->get_E(), 210e9);
    EXPECT_DOUBLE_EQ(elem->get_nu(), 0.3);
    EXPECT_DOUBLE_EQ(elem->get_rho(), 7850.0);
}

TEST_F(DataTest, ElementHasCorrectCoordinates) {

    Data data(parser2D);
    
    auto elem = data.get_elem(0);
    
    EXPECT_DOUBLE_EQ(elem->get_coord(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(elem->get_coord(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(elem->get_coord(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(elem->get_coord(1, 1), 0.0);
}

TEST_F(DataTest, NodesAreNumberedSequentially) {

    Data data(parser2D);
    
    for (int i = 0; i < data.nodes_count(); i++) {
        auto node = data.get_node(i);

        EXPECT_EQ(node->getID(), i + 1);
    }
}

class MockParserMultipleElements : public Parser {
public:
    MockParserMultipleElements() {
        settings.dimensions = "2D";
        settings.plane_state = "p-stress";
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
            0.5, 1.0, 0.0,
            1.0, 0.0, 0.0,
            2.0, 0.0, 0.0,
            1.5, 1.0, 0.0
        };

        mesh.elems_count = 2;
        mesh.elem_id = { 1, 2 };
        mesh.elem_types = { TRI, TRI };
        mesh.elem_orders = { 1, 1 };
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

TEST(DataMultipleElementsTest, CreatesMultipleElements) {

    auto parser = std::make_shared<MockParserMultipleElements>();
    Data data(parser);
    
    EXPECT_EQ(data.elements_count(), 2);
    EXPECT_EQ(data.get_elem(0)->get_type(), ElemType::TRI);
    EXPECT_EQ(data.get_elem(1)->get_type(), ElemType::TRI);
}

class MockParserSpectralQuad : public Parser {
public:
    MockParserSpectralQuad() {
        settings.dimensions = "2D";
        settings.plane_state = "p-stress";
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
            0.5, 0.0, 0.0,
            1.0, 0.5, 0.0,
            0.5, 1.0, 0.0,
            0.0, 0.5, 0.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { QUADSEM };
        mesh.elem_orders = { 3 };
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

TEST(DataSpectralQuadTest, SpectralQuad_CreatesCorrectNodeCount) {

    auto parser = std::make_shared<MockParserSpectralQuad>();
    Data data(parser);
    
    auto elem = data.get_elem(0);
    EXPECT_EQ(elem->get_order(), 3);
    EXPECT_EQ(elem->nodes_count(), 16);
}

TEST(DataSpectralQuadTest, SpectralQuad_CreatesGLLNodes) {

    auto parser = std::make_shared<MockParserSpectralQuad>();
    Data data(parser);
    
    EXPECT_GT(data.nodes_count(), 8);
}

TEST(DataSpectralQuadTest, SpectralQuad_ElementTypeIsCorrect) {

    auto parser = std::make_shared<MockParserSpectralQuad>();
    Data data(parser);
    
    EXPECT_EQ(data.get_elem(0)->get_type(), ElemType::QUADSEM);
}

class MockParserSpectralQuadConstrained : public Parser {
public:
    MockParserSpectralQuadConstrained() {
        settings.dimensions = "2D";
        settings.plane_state = "p-stress";
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
            0.5, 0.0, 0.0,
            1.0, 0.5, 0.0,
            0.5, 1.0, 0.0,
            0.0, 0.5, 0.0
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { QUADSEM };
        mesh.elem_orders = { 3 };
        mesh.elem_nodes = { 1, 2, 3, 4, 5, 6, 7, 8 };

        ::material mat;
        mat.id = 1;
        mat.type = 0;
        mat.constants = { 210e9, 0.3, 7850.0 };
        this->material.push_back(mat);

        ::restraints rest;
        rest.id = 1;
        rest.apply_to = { 1, 4, 8 };
        rest.size = 3;
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

TEST(DataSpectralQuadTest, SpectralQuadConstrained_PropagatesConstraints) {

    auto parser = std::make_shared<MockParserSpectralQuadConstrained>();
    Data data(parser);
    
    int constrained_count = 0;
    for (int i = 0; i < data.nodes_count(); i++) {
        if (!data.get_node(i)->constraints.empty()) {
            constrained_count++;
        }
    }
    
    EXPECT_GE(constrained_count, 3);
}

class MockParserSpectralHex : public Parser {
public:
    MockParserSpectralHex() {
        settings.dimensions = "3D";
        settings.plane_state = "";
        settings.analisys_type = "static";
        settings.d = 0.0;
        settings.max_time = 1.0;
        settings.max_iter = 100;
        settings.iter_res_output = 1;

        mesh.nodes_count = 20;
        mesh.node_id.resize(20);
        for (int i = 0; i < 20; i++) mesh.node_id[i] = i + 1;
        
        mesh.nodes_coord = {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            0.0, 1.0, 1.0,

            0.5, 0.0, 0.0,
            1.0, 0.5, 0.0,
            0.5, 1.0, 0.0,
            0.0, 0.5, 0.0,

            0.5, 0.0, 1.0,
            1.0, 0.5, 1.0,
            0.5, 1.0, 1.0,
            0.0, 0.5, 1.0,

            0.0, 0.0, 0.5,
            1.0, 0.0, 0.5,
            1.0, 1.0, 0.5,
            0.0, 1.0, 0.5
        };

        mesh.elems_count = 1;
        mesh.elem_id = { 1 };
        mesh.elem_types = { HEXSEM };
        mesh.elem_orders = { 3 };
        
        mesh.elem_nodes.resize(20);
        for (int i = 0; i < 20; i++) mesh.elem_nodes[i] = i + 1;

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

TEST(DataSpectralHexTest, SpectralHex_CreatesCorrectNodeCount) {

    auto parser = std::make_shared<MockParserSpectralHex>();
    Data data(parser);
    
    auto elem = data.get_elem(0);
    EXPECT_EQ(elem->get_order(), 3);
    EXPECT_EQ(elem->nodes_count(), 64);
}

TEST(DataSpectralHexTest, SpectralHex_CreatesGLLNodes) {

    auto parser = std::make_shared<MockParserSpectralHex>();
    Data data(parser);
    
    EXPECT_GT(data.nodes_count(), 20);
}

TEST(DataSpectralHexTest, SpectralHex_ElementTypeIsCorrect) {

    auto parser = std::make_shared<MockParserSpectralHex>();
    Data data(parser);
    
    EXPECT_EQ(data.get_elem(0)->get_type(), ElemType::HEXSEM);
}

class MockParserWithInfinite2D : public Parser {
public:
    MockParserWithInfinite2D() {
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

        ::infinite inf;
        inf.apply_to = { 1, 2 };
        inf.size = 1;
        inf.point = { 0.0, 0.5, 0.0 };
        this->Parser::infinite.push_back(inf);

        this->Parser::load.clear();
    }
};

TEST(DataInfiniteElementTest, InfiniteElement2D_CreatesInfQuad) {

    auto parser = std::make_shared<MockParserWithInfinite2D>();
    Data data(parser);
    
    EXPECT_EQ(data.elements_count(), 2);
}

TEST(DataInfiniteElementTest, InfiniteElement2D_CreatesAdditionalNodes) {

    auto parser = std::make_shared<MockParserWithInfinite2D>();
    Data data(parser);
    
    EXPECT_GT(data.nodes_count(), 4);
}

TEST(DataInfiniteElementTest, InfiniteElement2D_HasCorrectType) {

    auto parser = std::make_shared<MockParserWithInfinite2D>();
    Data data(parser);
    
    EXPECT_EQ(data.get_elem(1)->get_type(), ElemType::INFQUAD);
}

TEST(DataInfiniteElementTest, InfiniteElement2D_CountsInfElements) {

    auto parser = std::make_shared<MockParserWithInfinite2D>();
    Data data(parser);
    
    EXPECT_EQ(data.num_inf_elems, 1);
}

TEST(DataRenumberingTest, NodesAreSequential) {

    auto parser = std::make_shared<MockParser>();
    Data data(parser);
    
    for (int i = 0; i < data.nodes_count(); i++) {
        EXPECT_EQ(data.get_node(i)->getID(), i + 1);
    }
}

TEST(DataRenumberingTest, ElementNodesAreValid) {

    auto parser = std::make_shared<MockParser>();
    Data data(parser);
    
    for (int e = 0; e < data.elements_count(); e++) {
        auto elem = data.get_elem(e);
        for (int n = 0; n < elem->nodes_count(); n++) {
            int node_id = elem->get_node(n);
            EXPECT_GE(node_id, 1);
            EXPECT_LE(node_id, data.nodes_count());
        }
    }
}

TEST(DataRenumberingTest, ConstraintsPreservedAfterRenumbering) {

    auto parser = std::make_shared<MockParserWithConstraints>();
    Data data(parser);
    
    int constrained_count = 0;
    for (int i = 0; i < data.nodes_count(); i++) {
        if (!data.get_node(i)->constraints.empty()) {
            constrained_count++;
        }
    }
    
    EXPECT_EQ(constrained_count, 2);
}

TEST(DataEdgeCasesTest, SpectralQuadOrder4) {

    class MockParserQuadOrder4 : public Parser {
    public:
        MockParserQuadOrder4() {
            settings.dimensions = "2D";
            settings.plane_state = "p-stress";
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
                0.5, 0.0, 0.0,
                1.0, 0.5, 0.0,
                0.5, 1.0, 0.0,
                0.0, 0.5, 0.0
            };

            mesh.elems_count = 1;
            mesh.elem_id = { 1 };
            mesh.elem_types = { QUADSEM };
            mesh.elem_orders = { 4 };
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
    
    auto parser = std::make_shared<MockParserQuadOrder4>();
    Data data(parser);
    
    EXPECT_EQ(data.get_elem(0)->nodes_count(), 25);
}

TEST(DataEdgeCasesTest, SpectralQuadOrder5) {

    class MockParserQuadOrder5 : public Parser {
    public:
        MockParserQuadOrder5() {
            settings.dimensions = "2D";
            settings.plane_state = "p-stress";
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
                0.5, 0.0, 0.0,
                1.0, 0.5, 0.0,
                0.5, 1.0, 0.0,
                0.0, 0.5, 0.0
            };

            mesh.elems_count = 1;
            mesh.elem_id = { 1 };
            mesh.elem_types = { QUADSEM };
            mesh.elem_orders = { 5 };
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
    
    auto parser = std::make_shared<MockParserQuadOrder5>();
    Data data(parser);
    
    EXPECT_EQ(data.get_elem(0)->nodes_count(), 36);
}

TEST_F(DataTest, DMatrix2D_PlaneStress_ExactValues) {

    Data data(parser2D);
    
    auto elem = data.get_elem(0);
    double E = 210e9;
    double nu = 0.3;
    double factor = E / (1 - nu * nu);
    
    EXPECT_NEAR(elem->D(0, 0), factor * 1.0, 1e-6);
    EXPECT_NEAR(elem->D(1, 1), factor * 1.0, 1e-6);
    EXPECT_NEAR(elem->D(2, 2), factor * (1 - nu) / 2.0, 1e-6);
    
    EXPECT_NEAR(elem->D(0, 1), factor * nu, 1e-6);
    EXPECT_NEAR(elem->D(1, 0), factor * nu, 1e-6);
    
    EXPECT_NEAR(elem->D(0, 2), 0.0, 1e-9);
    EXPECT_NEAR(elem->D(1, 2), 0.0, 1e-9);
    EXPECT_NEAR(elem->D(2, 0), 0.0, 1e-9);
    EXPECT_NEAR(elem->D(2, 1), 0.0, 1e-9);
}

class MockParserPlaneStrain : public Parser {
public:
    MockParserPlaneStrain() {
        settings.dimensions = "2D";
        settings.plane_state = "p-strain";
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

TEST(DataDMatrixTest, DMatrix2D_PlaneStrain_ExactValues) {

    auto parser = std::make_shared<MockParserPlaneStrain>();
    Data data(parser);
    
    auto elem = data.get_elem(0);
    double E = 210e9;
    double nu = 0.3;
    double factor = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    
    EXPECT_NEAR(elem->D(0, 0), factor * 1.0, 1e-6);
    EXPECT_NEAR(elem->D(1, 1), factor * 1.0, 1e-6);
    EXPECT_NEAR(elem->D(2, 2), factor * (1 - 2 * nu) / (2 * (1 - nu)), 1e-6);
    
    double nu_eff = nu / (1 - nu);
    EXPECT_NEAR(elem->D(0, 1), factor * nu_eff, 1e-6);
    EXPECT_NEAR(elem->D(1, 0), factor * nu_eff, 1e-6);
}

TEST_F(DataTest, DMatrix3D_ExactValues) {

    Data data(parser3D);
    
    auto elem = data.get_elem(0);
    double E = 210e9;
    double nu = 0.3;
    double factor = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    
    EXPECT_NEAR(elem->D(0, 0), factor * 1.0, 1e-6);
    EXPECT_NEAR(elem->D(1, 1), factor * 1.0, 1e-6);
    EXPECT_NEAR(elem->D(2, 2), factor * 1.0, 1e-6);
    
    double shear_factor = factor * (1 - 2 * nu) / (2 * (1 - nu));
    EXPECT_NEAR(elem->D(3, 3), shear_factor, 1e-6);
    EXPECT_NEAR(elem->D(4, 4), shear_factor, 1e-6);
    EXPECT_NEAR(elem->D(5, 5), shear_factor, 1e-6);
    
    double nu_eff = nu / (1 - nu);
    EXPECT_NEAR(elem->D(0, 1), factor * nu_eff, 1e-6);
    EXPECT_NEAR(elem->D(0, 2), factor * nu_eff, 1e-6);
    EXPECT_NEAR(elem->D(1, 2), factor * nu_eff, 1e-6);
    
    EXPECT_NEAR(elem->D(1, 0), elem->D(0, 1), 1e-9);
    EXPECT_NEAR(elem->D(2, 0), elem->D(0, 2), 1e-9);
    EXPECT_NEAR(elem->D(2, 1), elem->D(1, 2), 1e-9);
}

TEST_F(DataTest, DMatrix_AllElementsHaveSameD) {

    class MockParserTwoElements : public Parser {
    public:
        MockParserTwoElements() {
            settings.dimensions = "2D";
            settings.plane_state = "p-stress";
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
                1.0, 1.0, 0.0,
                0.0, 1.0, 0.0,
                2.0, 0.0, 0.0,
                2.0, 1.0, 0.0
            };

            mesh.elems_count = 2;
            mesh.elem_id = { 1, 2 };
            mesh.elem_types = { QUAD, QUAD };
            mesh.elem_orders = { 1, 1 };
            mesh.elem_nodes = { 1, 2, 3, 4, 2, 5, 6, 3 };

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
    
    auto parser = std::make_shared<MockParserTwoElements>();
    Data data(parser);
    
    auto elem1 = data.get_elem(0);
    auto elem2 = data.get_elem(1);
    
    EXPECT_EQ(elem1->D.rows(), elem2->D.rows());
    EXPECT_EQ(elem1->D.cols(), elem2->D.cols());
    for (int i = 0; i < elem1->D.rows(); i++) {
        for (int j = 0; j < elem1->D.cols(); j++) {
            EXPECT_NEAR(elem1->D(i, j), elem2->D(i, j), 1e-9);
        }
    }
}

TEST(DataInfiniteElementTest, InfiniteElement2D_PoleCoordinates) {

    auto parser = std::make_shared<MockParserWithInfinite2D>();
    Data data(parser);
    
    auto inf_elem = data.get_elem(1);
    EXPECT_EQ(inf_elem->get_type(), ElemType::INFQUAD);
    
    auto* inf_quad = dynamic_cast<InfQuad*>(inf_elem.get());
    ASSERT_NE(inf_quad, nullptr);
    
    EXPECT_DOUBLE_EQ(inf_quad->pole_x, 0.0);
    EXPECT_DOUBLE_EQ(inf_quad->pole_y, 0.5);
    EXPECT_DOUBLE_EQ(inf_quad->pole_z, 0.0);
}

TEST(DataInfiniteElementTest, InfiniteElement2D_NodeCoordinates) {

    auto parser = std::make_shared<MockParserWithInfinite2D>();
    Data data(parser);
    
    auto inf_elem = data.get_elem(1);
    
    EXPECT_EQ(inf_elem->nodes_count(), 4);
    
    auto parent_elem = data.get_elem(0);
    std::vector<int> boundary_nodes = parent_elem->edge_to_node(2);
    
    int boundary_node_id = parent_elem->get_node(boundary_nodes[0]);
    auto boundary_node = data.get_node(boundary_node_id - 1);
    
    double x0 = inf_elem->get_coord(0, 0);
    double y0 = inf_elem->get_coord(0, 1);
    
    EXPECT_NEAR(x0, 1.0, 1e-6);
}

class MockParserWithInfinite3D : public Parser {
public:
    MockParserWithInfinite3D() {
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

        ::infinite inf;
        inf.apply_to = { 1, 1 };
        inf.size = 1;
        inf.point = { 2.0, 0.5, 0.5 };
        this->Parser::infinite.push_back(inf);

        this->restraints.clear();
        this->Parser::load.clear();
    }
};

TEST(DataInfiniteElementTest, InfiniteElement3D_CreatesInfHex) {

    auto parser = std::make_shared<MockParserWithInfinite3D>();
    Data data(parser);
    
    EXPECT_EQ(data.elements_count(), 2);
    EXPECT_EQ(data.get_elem(1)->get_type(), ElemType::INFHEX);
    EXPECT_EQ(data.num_inf_elems, 1);
}

TEST(DataInfiniteElementTest, InfiniteElement3D_PoleCoordinates) {

    auto parser = std::make_shared<MockParserWithInfinite3D>();
    Data data(parser);
    
    auto inf_elem = data.get_elem(1);
    auto* inf_hex = dynamic_cast<InfHex*>(inf_elem.get());
    ASSERT_NE(inf_hex, nullptr);
    
    EXPECT_DOUBLE_EQ(inf_hex->pole_x, 2.0);
    EXPECT_DOUBLE_EQ(inf_hex->pole_y, 0.5);
    EXPECT_DOUBLE_EQ(inf_hex->pole_z, 0.5);
}

TEST(DataSpectralQuadTest, SpectralQuad_GLLNodesCoordinates) {

    auto parser = std::make_shared<MockParserSpectralQuad>();
    Data data(parser);
    
    auto elem = data.get_elem(0);
    
    EXPECT_EQ(elem->nodes_count(), 16);
    
    for (int i = 0; i < elem->nodes_count(); i++) {
        double x = elem->get_coord(i, 0);
        double y = elem->get_coord(i, 1);
        
        EXPECT_GE(x, -0.1);
        EXPECT_LE(x, 1.1);
        EXPECT_GE(y, -0.1);
        EXPECT_LE(y, 1.1);
    }
}

TEST(DataSpectralQuadTest, SpectralQuad_OriginalNodesPreserved) {

    auto parser = std::make_shared<MockParserSpectralQuad>();
    Data data(parser);
    
    auto elem = data.get_elem(0);
    
    bool found_corner_0_0 = false;
    bool found_corner_1_1 = false;
    
    for (int i = 0; i < elem->nodes_count(); i++) {
        double x = elem->get_coord(i, 0);
        double y = elem->get_coord(i, 1);
        
        if (std::abs(x - 0.0) < 1e-6 && std::abs(y - 0.0) < 1e-6) {
            found_corner_0_0 = true;
        }
        if (std::abs(x - 1.0) < 1e-6 && std::abs(y - 1.0) < 1e-6) {
            found_corner_1_1 = true;
        }
    }
    
    EXPECT_TRUE(found_corner_0_0);
    EXPECT_TRUE(found_corner_1_1);
}

TEST(DataSpectralHexTest, SpectralHex_GLLNodesCoordinates) {

    auto parser = std::make_shared<MockParserSpectralHex>();
    Data data(parser);
    
    auto elem = data.get_elem(0);
    
    EXPECT_EQ(elem->nodes_count(), 64);
    
    for (int i = 0; i < elem->nodes_count(); i++) {
        double x = elem->get_coord(i, 0);
        double y = elem->get_coord(i, 1);
        double z = elem->get_coord(i, 2);
        
        EXPECT_GE(x, -0.1);
        EXPECT_LE(x, 1.1);
        EXPECT_GE(y, -0.1);
        EXPECT_LE(y, 1.1);
        EXPECT_GE(z, -0.1);
        EXPECT_LE(z, 1.1);
    }
}

TEST(DataRenumberingTest, Renumbering_NoGapsInSequence) {

    auto parser = std::make_shared<MockParser>();
    Data data(parser);
    
    std::set<int> node_ids;
    for (int i = 0; i < data.nodes_count(); i++) {
        node_ids.insert(data.get_node(i)->getID());
    }
    
    for (int i = 1; i <= data.nodes_count(); i++) {
        EXPECT_TRUE(node_ids.count(i) > 0) << "Missing node ID: " << i;
    }
}

TEST(DataRenumberingTest, Renumbering_ElementNodeReferencesValid) {

    class MockParserNonSequential : public Parser {
    public:
        MockParserNonSequential() {
            settings.dimensions = "2D";
            settings.plane_state = "p-stress";
            settings.analisys_type = "static";
            settings.d = 0.0;
            settings.max_time = 1.0;
            settings.max_iter = 100;
            settings.iter_res_output = 1;

            mesh.nodes_count = 4;
            mesh.node_id = { 10, 20, 30, 40 };
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
            mesh.elem_nodes = { 10, 20, 30, 40 };

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
    
    auto parser = std::make_shared<MockParserNonSequential>();
    Data data(parser);
    
    auto elem = data.get_elem(0);
    for (int i = 0; i < elem->nodes_count(); i++) {
        int node_id = elem->get_node(i);
        EXPECT_GE(node_id, 1);
        EXPECT_LE(node_id, data.nodes_count());
        
        auto node = data.get_node(node_id - 1);
        ASSERT_NE(node, nullptr);
        EXPECT_EQ(node->getID(), node_id);
    }
}

TEST(DataRenumberingTest, Renumbering_ConstraintsPreserved) {

    auto parser = std::make_shared<MockParserWithConstraints>();
    Data data(parser);
    
    int constrained_count = 0;
    std::map<int, std::map<int, double>> constraints_map;
    
    for (int i = 0; i < data.nodes_count(); i++) {
        auto node = data.get_node(i);
        if (!node->constraints.empty()) {
            constrained_count++;
            constraints_map[node->getID()] = node->constraints;
        }
    }
    
    EXPECT_EQ(constrained_count, 2);
    
    for (const auto& [node_id, constraints] : constraints_map) {
        EXPECT_EQ(constraints.count(0), 1);
        EXPECT_EQ(constraints.count(1), 1);
        EXPECT_DOUBLE_EQ(constraints.at(0), 0.0);
        EXPECT_DOUBLE_EQ(constraints.at(1), 0.0);
    }
}

TEST(DataRenumberingTest, Renumbering_NodeCoordinatesPreserved) {

    auto parser = std::make_shared<MockParser>();
    Data data(parser);
    
    EXPECT_DOUBLE_EQ(data.get_node(0)->getX(), 0.0);
    EXPECT_DOUBLE_EQ(data.get_node(0)->getY(), 0.0);
    EXPECT_DOUBLE_EQ(data.get_node(1)->getX(), 1.0);
    EXPECT_DOUBLE_EQ(data.get_node(1)->getY(), 0.0);
    EXPECT_DOUBLE_EQ(data.get_node(2)->getX(), 1.0);
    EXPECT_DOUBLE_EQ(data.get_node(2)->getY(), 1.0);
    EXPECT_DOUBLE_EQ(data.get_node(3)->getX(), 0.0);
    EXPECT_DOUBLE_EQ(data.get_node(3)->getY(), 1.0);
}

TEST_F(DataConstraintsTest, Constraints_AppliedToCorrectNodes) {

    Data data(parser);
    
    std::vector<int> constrained_node_ids;
    for (int i = 0; i < data.nodes_count(); i++) {
        auto node = data.get_node(i);
        if (!node->constraints.empty()) {
            constrained_node_ids.push_back(node->getID());
        }
    }
    
    EXPECT_EQ(constrained_node_ids.size(), 2);
}

TEST_F(DataConstraintsTest, Constraints_CorrectDOFValues) {

    Data data(parser);
    
    for (int i = 0; i < data.nodes_count(); i++) {
        auto node = data.get_node(i);
        if (!node->constraints.empty()) {

            EXPECT_TRUE(node->constraints.count(0) > 0);
            EXPECT_TRUE(node->constraints.count(1) > 0);
            
            EXPECT_DOUBLE_EQ(node->constraints.at(0), 0.0);
            EXPECT_DOUBLE_EQ(node->constraints.at(1), 0.0);
            
            EXPECT_TRUE(node->constraints.count(2) == 0);
        }
    }
}

class MockParserWithConstraints3D : public Parser {
public:
    MockParserWithConstraints3D() {
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

        ::restraints rest;
        rest.id = 1;
        rest.apply_to = { 1, 2, 3, 4 };
        rest.size = 4;
        rest.cs = 0;
        rest.data[0] = 0.0;
        rest.data[1] = 0.0;
        rest.data[2] = 0.0;
        rest.flag[0] = 1;
        rest.flag[1] = 1;
        rest.flag[2] = 1;
        rest.flag[3] = 0;
        rest.flag[4] = 0;
        rest.flag[5] = 0;
        this->restraints.push_back(rest);

        this->Parser::load.clear();
        this->Parser::infinite.clear();
    }
};

TEST(DataConstraintsTest, Constraints3D_AppliedCorrectly) {

    auto parser = std::make_shared<MockParserWithConstraints3D>();
    Data data(parser);
    
    int constrained_count = 0;
    for (int i = 0; i < data.nodes_count(); i++) {
        auto node = data.get_node(i);
        if (!node->constraints.empty()) {
            constrained_count++;

            EXPECT_TRUE(node->constraints.count(0) > 0);
            EXPECT_TRUE(node->constraints.count(1) > 0);
            EXPECT_TRUE(node->constraints.count(2) > 0);
        }
    }
    
    EXPECT_EQ(constrained_count, 4);
}

TEST_F(DataConstraintsTest, Constraints_OnlyFlaggedDOFsApplied) {

    Data data(parser);
    
    for (int i = 0; i < data.nodes_count(); i++) {
        auto node = data.get_node(i);
        if (!node->constraints.empty()) {

            EXPECT_TRUE(node->constraints.count(2) == 0);
        }
    }
}

#endif
