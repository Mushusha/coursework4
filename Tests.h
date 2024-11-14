#pragma once
#include "Parser.h"
#include "Data.h"
#include "Element.h"
#include "TriangleElements.h"


#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"


void test_triElement_B();
void test_triElement_locK();
void test_triElement_globK();
void test_triElement_locR();
void test_triElement_globR();
void test_triElement_constraintsK();
void test_triElement_constraintsR();


void tests();