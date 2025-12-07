#pragma once

#include <gtest/gtest.h>
#include <filesystem>
#include <cmath>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "Calc.h"

namespace fs = std::filesystem;

class IntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        tolerance = 1e-9;
    }

    double tolerance;

    bool compareVTUFiles(const std::string& resultFile, const std::string& expectedFile) {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> resultReader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        resultReader->SetFileName(resultFile.c_str());
        resultReader->Update();

        vtkSmartPointer<vtkXMLUnstructuredGridReader> expectedReader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        expectedReader->SetFileName(expectedFile.c_str());
        expectedReader->Update();

        vtkUnstructuredGrid* resultGrid = resultReader->GetOutput();
        vtkUnstructuredGrid* expectedGrid = expectedReader->GetOutput();

        if (!resultGrid || !expectedGrid) {
            std::cerr << "Failed to read VTU files" << std::endl;
            return false;
        }

        if (resultGrid->GetNumberOfPoints() != expectedGrid->GetNumberOfPoints()) {
            std::cerr << "Point count mismatch: " << resultGrid->GetNumberOfPoints()
                << " vs " << expectedGrid->GetNumberOfPoints() << std::endl;
            return false;
        }

        if (resultGrid->GetNumberOfCells() != expectedGrid->GetNumberOfCells()) {
            std::cerr << "Cell count mismatch: " << resultGrid->GetNumberOfCells()
                << " vs " << expectedGrid->GetNumberOfCells() << std::endl;
            return false;
        }

        if (!comparePoints(resultGrid, expectedGrid)) {
            return false;
        }

        if (!comparePointData(resultGrid, expectedGrid)) {
            return false;
        }

        if (!compareCellData(resultGrid, expectedGrid)) {
            return false;
        }

        return true;
    }

    bool comparePoints(vtkUnstructuredGrid* result, vtkUnstructuredGrid* expected) {
        vtkIdType numPoints = result->GetNumberOfPoints();
        for (vtkIdType i = 0; i < numPoints; ++i) {
            double resultPt[3], expectedPt[3];
            result->GetPoint(i, resultPt);
            expected->GetPoint(i, expectedPt);

            for (int j = 0; j < 3; ++j) {
                if (std::abs(resultPt[j] - expectedPt[j]) > tolerance) {
                    std::cerr << "Point " << i << " coordinate " << j << " mismatch: "
                        << resultPt[j] << " vs " << expectedPt[j] << std::endl;
                    return false;
                }
            }
        }
        return true;
    }

    bool comparePointData(vtkUnstructuredGrid* result, vtkUnstructuredGrid* expected) {
        vtkPointData* resultData = result->GetPointData();
        vtkPointData* expectedData = expected->GetPointData();

        if (resultData->GetNumberOfArrays() != expectedData->GetNumberOfArrays()) {
            std::cerr << "Point data array count mismatch" << std::endl;
            return false;
        }

        for (int a = 0; a < expectedData->GetNumberOfArrays(); ++a) {
            vtkDataArray* expectedArray = expectedData->GetArray(a);
            vtkDataArray* resultArray = resultData->GetArray(expectedArray->GetName());

            if (!resultArray) {
                std::cerr << "Missing point data array: " << expectedArray->GetName() << std::endl;
                return false;
            }

            if (!compareDataArray(resultArray, expectedArray, expectedArray->GetName())) {
                return false;
            }
        }

        return true;
    }

    bool compareCellData(vtkUnstructuredGrid* result, vtkUnstructuredGrid* expected) {
        vtkCellData* resultData = result->GetCellData();
        vtkCellData* expectedData = expected->GetCellData();

        if (resultData->GetNumberOfArrays() != expectedData->GetNumberOfArrays()) {
            std::cerr << "Cell data array count mismatch" << std::endl;
            return false;
        }

        for (int a = 0; a < expectedData->GetNumberOfArrays(); ++a) {
            vtkDataArray* expectedArray = expectedData->GetArray(a);
            vtkDataArray* resultArray = resultData->GetArray(expectedArray->GetName());

            if (!resultArray) {
                std::cerr << "Missing cell data array: " << expectedArray->GetName() << std::endl;
                return false;
            }

            if (!compareDataArray(resultArray, expectedArray, expectedArray->GetName())) {
                return false;
            }
        }

        return true;
    }

    bool compareDataArray(vtkDataArray* result, vtkDataArray* expected, const char* name) {
        if (result->GetNumberOfTuples() != expected->GetNumberOfTuples()) {
            std::cerr << "Array " << name << " tuple count mismatch: "
                << result->GetNumberOfTuples() << " vs " << expected->GetNumberOfTuples() << std::endl;
            return false;
        }

        if (result->GetNumberOfComponents() != expected->GetNumberOfComponents()) {
            std::cerr << "Array " << name << " component count mismatch: "
                << result->GetNumberOfComponents() << " vs " << expected->GetNumberOfComponents() << std::endl;
            return false;
        }

        vtkIdType numTuples = expected->GetNumberOfTuples();
        int numComponents = expected->GetNumberOfComponents();

        for (vtkIdType i = 0; i < numTuples; ++i) {
            for (int c = 0; c < numComponents; ++c) {
                double resultVal = result->GetComponent(i, c);
                double expectedVal = expected->GetComponent(i, c);

                double maxVal = std::max(std::abs(resultVal), std::abs(expectedVal));
                double effectiveTolerance = (maxVal > 1.0) ? tolerance * maxVal : tolerance;

                if (std::abs(resultVal - expectedVal) > effectiveTolerance) {
                    std::cerr << "Array " << name << " value mismatch at tuple " << i
                        << " component " << c << ": "
                        << resultVal << " vs " << expectedVal << std::endl;
                    return false;
                }
            }
        }

        return true;
    }
};

class Kq18000IntegrationTest : public IntegrationTest {
protected:
    void SetUp() override {
        IntegrationTest::SetUp();

        testDataDir = "../Tests/data/";
        
        inputFile = testDataDir + "kq18000";
        expectedOutputFile = testDataDir + "kq18000_expected.vtu";
        actualOutputFile = testDataDir + "kq18000.vtu";
    }

    void TearDown() override {

        if (fs::exists(actualOutputFile)) {
            fs::remove(actualOutputFile);
        }
    }

    std::string testDataDir;
    std::string inputFile;
    std::string expectedOutputFile;
    std::string actualOutputFile;
};

