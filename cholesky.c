#include "cholesky.h"
#include "matrixutils.h"

#include <stdio.h>
#include <math.h>

bool BuildCholeskyLT (double **A, int matrixSize, int *D)
{
    for (int step = 0; step < matrixSize; ++step)
    {
        double modifier;
        if (A[step][step] < 0)
        {
            D[step] = -1;
            modifier = -1.0 / sqrt (-A[step][step]);
        }
        else
        {
            modifier = 1.0 / sqrt (A[step][step]);
            D[step] = 1;
        }

        MultiplyRowPart (A, matrixSize, matrixSize, step, modifier, step);
        for (int row = step + 1; row < matrixSize; ++row)
        {
            modifier = -A[row][step] / A[step][step];
            AddMultipliedRowPart (A, matrixSize, matrixSize, row, step, modifier, step);
        }
    }

    return true;
}

bool SolveCholesky (double **LT, int matrixSize, int *D, double **B, int results, double ***X)
{
    // Calculate L * Y = B.
    *X = CopyMatrix (B, matrixSize, results);
    for (int step = 0; step < matrixSize; ++step)
    {
        MultiplyRow (*X, matrixSize, results, step, 1.0 / LT[step][step]);
        for (int row = step + 1; row < matrixSize; ++row)
        {
            AddMultipliedRow (*X, matrixSize, results, row, step, -LT[step][row]);
        }
    }

    // Calculate D * LT * X = Y.
    for (int step = matrixSize - 1; step >= 0;--step)
    {
        MultiplyRow (*X, matrixSize, results, step, 1.0 / (D[step] * LT[step][step]));
        for (int row = step - 1; row >= 0; --row)
        {
            AddMultipliedRow (*X, matrixSize, results, row, step, -(D[row] * LT[row][step]));
        }
    }

    return true;
}
