#pragma once
#include <stdbool.h>

bool BuildCholeskyLT (double **A, int matrixSize, int *D);
bool SolveCholesky (double **LT, int matrixSize, int *D, double **B, int results, double ***X);
