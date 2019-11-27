#pragma once
#include <stdbool.h>

/// Note: B should be [matrixSize]x1 sized matrix.
bool SolveGMRESArnoldi (double **A, int matrixSize, double **B, double ***X);
