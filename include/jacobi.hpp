#pragma once
#include <vector>

using matrixType = std::vector <std::vector <double>>;

std::vector <double> jacobi(matrixType A, std::vector <double> b, std::vector <double> x, int dim);
void printMatrix(matrixType matrix,int N);
void printVector(std::vector <double> vector);