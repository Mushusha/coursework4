// #define ENABLE_INTEGRATION_TESTS

#include "StaticQuadTests.h"

#ifdef ENABLE_INTEGRATION_TESTS

TEST_F(Kq18000IntegrationTest, CalculationProducesCorrectResults) {

    std::string inputFileWithExt = inputFile + ".fc";
    ASSERT_TRUE(fs::exists(inputFileWithExt)) << "Input file " << inputFileWithExt << " does not exist";

    ASSERT_TRUE(fs::exists(expectedOutputFile)) 
        << "Expected output file " << expectedOutputFile << " does not exist";

    Calculate solve(inputFile);
    ASSERT_NO_THROW(solve.Solve()) << "Calculation threw an exception";

    ASSERT_TRUE(fs::exists(actualOutputFile)) 
        << "Output file " << actualOutputFile << " was not created";

    EXPECT_TRUE(compareVTUFiles(actualOutputFile, expectedOutputFile))
        << "VTU file comparison failed";
}

TEST_F(Kq18000IntegrationTest, OutputFileHasExpectedStructure) {

    if (!fs::exists(expectedOutputFile)) {
        GTEST_SKIP() << "Expected output file not found, skipping structure test";
    }

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(expectedOutputFile.c_str());
    reader->Update();

    vtkUnstructuredGrid* grid = reader->GetOutput();
    ASSERT_NE(grid, nullptr) << "Failed to read expected VTU file";

    EXPECT_EQ(grid->GetNumberOfPoints(), 18951) << "Unexpected number of points";
    EXPECT_EQ(grid->GetNumberOfCells(), 18671) << "Unexpected number of cells";

    vtkPointData* pointData = grid->GetPointData();
    ASSERT_NE(pointData, nullptr);

    EXPECT_NE(pointData->GetArray("Displacement"), nullptr) 
        << "Missing Displacement array";
    EXPECT_NE(pointData->GetArray("Stress"), nullptr) 
        << "Missing Stress array";
    EXPECT_NE(pointData->GetArray("nodeID"), nullptr) 
        << "Missing nodeID array";

    vtkDataArray* dispArray = pointData->GetArray("Displacement");
    if (dispArray) {
        EXPECT_EQ(dispArray->GetNumberOfComponents(), 2) 
            << "Displacement should have 2 components";
    }

    vtkDataArray* stressArray = pointData->GetArray("Stress");
    if (stressArray) {
        EXPECT_EQ(stressArray->GetNumberOfComponents(), 3) 
            << "Stress should have 3 components";
    }

    vtkCellData* cellData = grid->GetCellData();
    ASSERT_NE(cellData, nullptr);

    EXPECT_NE(cellData->GetArray("ElementID"), nullptr) 
        << "Missing ElementID array";
    EXPECT_NE(cellData->GetArray("ElementType"), nullptr) 
        << "Missing ElementType array";
}

TEST_F(Kq18000IntegrationTest, DisplacementValuesInExpectedRange) {
    if (!fs::exists(expectedOutputFile)) {
        GTEST_SKIP() << "Expected output file not found";
    }

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(expectedOutputFile.c_str());
    reader->Update();

    vtkUnstructuredGrid* grid = reader->GetOutput();
    vtkDataArray* dispArray = grid->GetPointData()->GetArray("Displacement");
    
    ASSERT_NE(dispArray, nullptr);

    double range[2];
    dispArray->GetRange(range, -1);

    EXPECT_GE(range[0], 0.0) << "Minimum displacement magnitude should be >= 0";
    EXPECT_LE(range[1], 1e-8) << "Maximum displacement is too large for this problem";
    
    EXPECT_NEAR(range[0], 4.537e-12, 1e-11) << "Minimum displacement out of expected range";
    EXPECT_NEAR(range[1], 4.949e-10, 1e-9) << "Maximum displacement out of expected range";
}

TEST_F(Kq18000IntegrationTest, StressValuesInExpectedRange) {
    if (!fs::exists(expectedOutputFile)) {
        GTEST_SKIP() << "Expected output file not found";
    }

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(expectedOutputFile.c_str());
    reader->Update();

    vtkUnstructuredGrid* grid = reader->GetOutput();
    vtkDataArray* stressArray = grid->GetPointData()->GetArray("Stress");
    
    ASSERT_NE(stressArray, nullptr);

    double range[2];
    stressArray->GetRange(range, -1);

    EXPECT_NEAR(range[0], 0.445, 0.01) << "Minimum stress out of expected range";
    EXPECT_NEAR(range[1], 29.265, 0.01) << "Maximum stress out of expected range";
}

#endif

