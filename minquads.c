#include "minquads.h"
#include "matrixutils.h"
#include "cholesky.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

bool SolveMinQuads (double **A, int rows, int cols, double **B, int results, double ***X)
{
    double **At = TransposeMatrix (A, rows, cols);
    double **AtxA = AllocateMatrix (cols, cols);
    double **AtxB = AllocateMatrix (cols, results);
    int *D = calloc (cols, sizeof (int));

    MultiplyMatrices (At, A, AtxA, cols, rows, cols);
    MultiplyMatrices (At, B, AtxB, cols, rows, results);
    BuildCholeskyLT (AtxA, cols, D);
    SolveCholesky (AtxA, cols, D, AtxB, results, X);

    FreeMatrix (At, cols, rows);
    FreeMatrix (AtxA, cols, cols);
    free (D);
    return true;
}
