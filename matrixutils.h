#pragma once

double **AllocateMatrix (int rows, int cols);
double **CopyMatrix (double **matrix, int rows, int cols);
void CopyMatrixInto (double **matrix, int rows, int cols, double **output);
double **TransformMatrixByRowOrder (double **matrix, int rows, int cols, int *rowOrder);
double **TransformMatrixByColOrder (double **matrix, int rows, int cols, int *colOrder);
double **TransposeMatrix (double **matrix, int rows, int cols);
void FreeMatrix (double **matrix, int rows, int cols);

double ColumnQuadricNormal (double **matrix, int rows, int colIndex);
double ColumnsScalarMultiplication (double **matrix, int rows, int col1Index, int col2Index);
void PrintMatrix (double **matrix, int rows, int cols);

void MultiplyRow (double **matrix, int rows, int cols, int row, double by);
void MultiplyRowPart (double **matrix, int rows, int cols, int row, double by, int start);
void MultiplyCol (double **matrix, int rows, int cols, int col, double by);
void SwapRows (double **matrix, int rows, int cols, int row1, int row2);
void SwapCols (double **matrix, int rows, int cols, int col1, int col2);

void MultiplyMatrices (double **first, double **second, double **output, int firstRows, int firstCols, int secondsCols);
void AddMultipliedRow (double **matrix, int rows, int cols, int dst, int src, double modifier);
void AddMultipliedRowPart (double **matrix, int rows, int cols, int dst, int src, double modifier, int startFrom);
void AddMultipliedCol (double **matrix, int rows, int cols, int dst, int src, double modifier);
