#include "VTUWriter.h"

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
        cell = vtkSmartPointer<vtkQuad>::New();
        break;
    case ElemType::TETRA:
        cell = vtkSmartPointer<vtkTetra>::New();
        break;
    case ElemType::HEX:
        cell = vtkSmartPointer<vtkHexahedron>::New();
        break;
    case ElemType::WEDGE:
        cell = vtkSmartPointer<vtkWedge>::New();
        break;
    case ElemType::PYR:
        cell = vtkSmartPointer<vtkPyramid>::New();
        break;
    case ElemType::QUADSEM: {
        int NODES = static_cast<int>(std::sqrt(nodeIds.size()));

        if (NODES * NODES != static_cast<int>(nodeIds.size()))
            throw std::runtime_error("VTUWriter: QUADSEM node count is not a perfect square");

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

    for (const auto& elem : elements)
        elemTypeArray->InsertNextValue(
            elem->get_type() == ElemType::QUADSEM
            ? static_cast<int>(ElemType::QUAD)
            : static_cast<int>(elem->get_type())
        );

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
