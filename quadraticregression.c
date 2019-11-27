#include "quadraticregression.h"
#include "matrixutils.h"
#include "householder.h"

#include <stdlib.h>
#include <math.h>

void QuadraticRegression (const double *x, const double *y, int count, int power, double ***output)
{
    double *values = calloc (power + power + 1, sizeof (double));
    for (int index = 0; index < power + power + 1; ++index)
    {
        double value = 0.0;
        for (int xIndex = 0; xIndex < count; ++xIndex)
        {
            value += pow (x[xIndex], index);
        }

        values[index] = value;
    }

    double **A = AllocateMatrix (power + 1, power + 1);
    *output = AllocateMatrix (power + 1, 1);

    for (int row = 0; row <= power; ++row)
    {
        for (int col = 0; col <= power; ++col)
        {
            A[row][col] = values[row + col];
        }
    }

    for (int row = 0; row <= power; ++row)
    {
        double value = 0.0;
        for (int index = 0; index < count; ++index)
        {
            value += y[index] * pow (x[index], row);
        }

        (*output)[row][0] = value;
    }

    SolveHouseholder (A, power + 1, *output, 1);
    FreeMatrix (A, power + 1, power + 1);
}
