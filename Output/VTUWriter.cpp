#include "VTUWriter.h"
#include <array>
#include <functional>

namespace {
    std::vector<int> buildHexSemToVtkMapping(int NODES) {
        const int N = NODES - 1;  // order
        const int nPts = NODES * NODES * NODES;
        std::vector<int> vtkToOurs(nPts);

        auto ourIdx = [NODES](int i, int j, int k) {
            return i + j * NODES + k * NODES * NODES;
        };

        int vtkIdx = 0;

        std::array<std::array<int, 3>, 8> corners = {{
            {0, 0, 0}, {N, 0, 0}, {N, N, 0}, {0, N, 0},
            {0, 0, N}, {N, 0, N}, {N, N, N}, {0, N, N}
        }};
        for (const auto& c : corners)
            vtkToOurs[vtkIdx++] = ourIdx(c[0], c[1], c[2]);

        std::array<std::pair<std::array<int, 3>, std::array<int, 3>>, 12> edges = {{
            {{0, 0, 0}, {N, 0, 0}}, {{N, 0, 0}, {N, N, 0}}, {{0, N, 0}, {N, N, 0}}, {{0, 0, 0}, {0, N, 0}},
            {{0, 0, N}, {N, 0, N}}, {{N, 0, N}, {N, N, N}}, {{0, N, N}, {N, N, N}}, {{0, 0, N}, {0, N, N}},
            {{0, 0, 0}, {0, 0, N}}, {{N, 0, 0}, {N, 0, N}}, {{N, N, 0}, {N, N, N}}, {{0, N, 0}, {0, N, N}}
        }};
        const int nEdge = (N > 1) ? (N - 1) : 0;
        for (const auto& e : edges) {
            for (int t = 1; t <= nEdge; ++t) {
                int i = e.first[0] + t * (e.second[0] - e.first[0]) / N;
                int j = e.first[1] + t * (e.second[1] - e.first[1]) / N;
                int k = e.first[2] + t * (e.second[2] - e.first[2]) / N;
                vtkToOurs[vtkIdx++] = ourIdx(i, j, k);
            }
        }

        std::array<std::function<void()>, 6> faceFns = {{
            [&]() { for (int v = 1; v <= nEdge; ++v) for (int u = 1; u <= nEdge; ++u) vtkToOurs[vtkIdx++] = ourIdx(0, u, v); },
            [&]() { for (int v = 1; v <= nEdge; ++v) for (int u = 1; u <= nEdge; ++u) vtkToOurs[vtkIdx++] = ourIdx(N, u, v); },
            [&]() { for (int v = 1; v <= nEdge; ++v) for (int u = 1; u <= nEdge; ++u) vtkToOurs[vtkIdx++] = ourIdx(u, 0, v); },
            [&]() { for (int v = 1; v <= nEdge; ++v) for (int u = 1; u <= nEdge; ++u) vtkToOurs[vtkIdx++] = ourIdx(u, N, v); },
            [&]() { for (int v = 1; v <= nEdge; ++v) for (int u = 1; u <= nEdge; ++u) vtkToOurs[vtkIdx++] = ourIdx(u, v, 0); },
            [&]() { for (int v = 1; v <= nEdge; ++v) for (int u = 1; u <= nEdge; ++u) vtkToOurs[vtkIdx++] = ourIdx(u, v, N); }
        }};
        for (const auto& fn : faceFns) fn();

        for (int k = 1; k <= nEdge; ++k) {
            for (int j = 1; j <= nEdge; ++j) {
                for (int i = 1; i <= nEdge; ++i) {
                    vtkToOurs[vtkIdx++] = ourIdx(i, j, k);
                }
            }
        }

        return vtkToOurs;
    }
}

bool VTUWriter::write(const std::string& filename) {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

    createPoints(unstructuredGrid);
    createCells(unstructuredGrid);
    addNodeData(unstructuredGrid);
    addElementData(unstructuredGrid);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writer->SetFileName(filename.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->SetDataModeToBinary();
    writer->SetCompressorTypeToZLib();

    logger& log = logger::log();
    log.print("End writing results");

    return writer->Write() == 1;
}
void VTUWriter::createPoints(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();

    for (const auto& node : nodes)
        points->InsertNextPoint(node->getX(), node->getY(), node->getZ());

    unstructuredGrid->SetPoints(points);
}

void VTUWriter::createCells(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
    for (const auto& elem : elements) {
        vtkSmartPointer<vtkCell> cell = createCell(elem->get_type(), elem->get_node());
        if (cell)
            unstructuredGrid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
    }
}

vtkSmartPointer<vtkCell> VTUWriter::createCell(ElemType type, const std::vector<int>& nodeIds) {
    vtkSmartPointer<vtkCell> cell;

    switch (type) {
    case ElemType::TRI:
        cell = vtkSmartPointer<vtkTriangle>::New();
        break;
    case ElemType::QUAD:
    case ElemType::INFQUAD:
        cell = vtkSmartPointer<vtkQuad>::New();
        break;
    case ElemType::TETRA:
        cell = vtkSmartPointer<vtkTetra>::New();
        break;
    case ElemType::HEX:
    case ElemType::INFHEX:
        cell = vtkSmartPointer<vtkHexahedron>::New();
        break;
    case ElemType::WEDGE:
        cell = vtkSmartPointer<vtkWedge>::New();
        break;
    case ElemType::PYR:
        cell = vtkSmartPointer<vtkPyramid>::New();
        break;
    case ElemType::QUADSEM:
    //case ElemType::INFQUADSEM: 
    {
        int NODES = static_cast<int>(std::sqrt(nodeIds.size()));

        if (NODES * NODES != static_cast<int>(nodeIds.size()))
            throw std::runtime_error("VTUWriter: QUADSEM/INFQUADSEM node count is not a perfect square");

        std::vector<int> boundary;

        for (int i = 0; i < NODES; ++i)
            boundary.push_back(nodeIds[i]);

        for (int j = 1; j < NODES - 1; ++j)
            boundary.push_back(nodeIds[j * NODES + (NODES - 1)]);

        for (int i = NODES - 1; i >= 0; --i)
            boundary.push_back(nodeIds[(NODES - 1) * NODES + i]);

        for (int j = NODES - 2; j >= 1; --j)
            boundary.push_back(nodeIds[j * NODES]);

        vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
        polygon->GetPointIds()->SetNumberOfIds(boundary.size());

        for (size_t i = 0; i < boundary.size(); ++i)
            polygon->GetPointIds()->SetId(i, boundary[i] - 1);

        return polygon;
    }
    case ElemType::HEXSEM:
    //case ElemType::INFHEXSEM: 
    {
        int NODES = static_cast<int>(std::cbrt(nodeIds.size()));
        int nPts = static_cast<int>(nodeIds.size());

        if (NODES * NODES * NODES != nPts)
            throw std::runtime_error("VTUWriter: HEXSEM/INFHEXSEM node count is not a perfect cube");

        std::vector<int> vtkToOurs = buildHexSemToVtkMapping(NODES);

        vtkSmartPointer<vtkLagrangeHexahedron> hex = vtkSmartPointer<vtkLagrangeHexahedron>::New();
        hex->SetUniformOrderFromNumPoints(nPts);
        hex->GetPointIds()->SetNumberOfIds(nPts);
        for (int vtkIdx = 0; vtkIdx < nPts; ++vtkIdx) {
            int ourIdx = vtkToOurs[vtkIdx];
            hex->GetPointIds()->SetId(vtkIdx, nodeIds[ourIdx] - 1);
        }

        return hex;
    }
    default:
        return nullptr;
    }

    for (size_t i = 0; i < nodeIds.size(); ++i)
        cell->GetPointIds()->SetId(i, nodeIds[i] - 1);

    return cell;
}

void VTUWriter::addNodeData(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
    if (nodes.empty() || nodes[0]->results.empty()) return;

    for (const auto& result_pair : nodes[0]->results) {
        const std::string& res_name = toString(result_pair.first);
        const size_t num_components = result_pair.second.size();

        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(res_name.c_str());
        dataArray->SetNumberOfComponents(num_components);

        for (const auto& node : nodes) {
            const auto& res_values = node->results.at(result_pair.first);
            dataArray->InsertNextTuple(res_values.data());
        }

        unstructuredGrid->GetPointData()->AddArray(dataArray);
    }

    vtkSmartPointer<vtkIntArray> nodeIDArray = vtkSmartPointer<vtkIntArray>::New();
    nodeIDArray->SetName("nodeID");
    nodeIDArray->SetNumberOfComponents(1);

    for (size_t i = 0; i < nodes.size(); ++i)
        nodeIDArray->InsertNextValue(i + 1);

    unstructuredGrid->GetPointData()->AddArray(nodeIDArray);
}

void VTUWriter::addElementData(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
    vtkSmartPointer<vtkIntArray> elemIdArray =
        vtkSmartPointer<vtkIntArray>::New();
    elemIdArray->SetName("ElementID");
    for (const auto& elem : elements)
        elemIdArray->InsertNextValue(elem->get_id());

    unstructuredGrid->GetCellData()->AddArray(elemIdArray);

    vtkSmartPointer<vtkIntArray> elemTypeArray =
        vtkSmartPointer<vtkIntArray>::New();
    elemTypeArray->SetName("ElementType");

    for (const auto& elem : elements) {
        int vtk_type;
        switch (elem->get_type()) {
        case ElemType::QUADSEM:
        case ElemType::INFQUADSEM:
            vtk_type = static_cast<int>(ElemType::QUAD);
            break;
        case ElemType::HEXSEM:
        case ElemType::INFHEXSEM:
            vtk_type = static_cast<int>(ElemType::HEX);
            break;
        default:
            vtk_type = static_cast<int>(elem->get_type());
            break;
        }
        elemTypeArray->InsertNextValue(vtk_type);
    }

    unstructuredGrid->GetCellData()->AddArray(elemTypeArray);
}

std::string VTUWriter::toString(ResType type) {
    switch (type) {
        case ResType::DISPLACEMENT: return "Displacement";
        case ResType::STRESS: return "Stress";
        case ResType::STRAIN: return "Strain";
        default: return "Unknown";
    }
}
