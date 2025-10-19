#pragma once

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkPyramid.h>
#include <vtkVertex.h>
#include <vtkPolygon.h>

#include "Elements/Element.h"
#include "Node/Node.h"

class VTUWriter {
public:
    VTUWriter(const std::vector<std::shared_ptr<Element>> elements,
                 const std::vector<std::shared_ptr<Node>> nodes)
        : elements(elements), nodes(nodes) {
    }

    bool write(const std::string& filename);

private:
    const std::vector<std::shared_ptr<Element>> elements;
    const std::vector<std::shared_ptr<Node>> nodes;

    void createPoints(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid);
    void createCells(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid);

    vtkSmartPointer<vtkCell> createCell(ElemType type, const std::vector<int>& nodeIds);

    void addNodeData(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid);
    void addElementData(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid);

    std::string toString(ResType type);
};