#include "householder.h"
#include "matrixutils.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double FindNormal (double **A, int matrixSize, int step)
{
    double normal = 0.0;
    for (int index = step; index < matrixSize; ++index)
    {
        normal += A[index][step] * A[index][step];
    }

    return sqrt (normal);
}

static double FindVectorNormal (double *v, int size)
{
    double normal = 0.0;
    for (int index = 0; index < size; ++index)
    {
        normal += v[index] * v[index];
    }

    return sqrt (normal);
}

static void ApplyW (double **matrix, int rows, int col, int step, double *w)
{
    double scalar = 0.0;
    for (int index = step; index < rows; ++index)
    {
        scalar += matrix[index][col] * w[index];
    }

    scalar *= 2.0;
    for (int index = step; index < rows; ++index)
    {
        matrix[index][col] -= scalar * w[index];
    }
}

bool SolveHouseholder (double **A, int matrixSize, double **B, int results)
{
    double *w = calloc (matrixSize, sizeof (double));
    for (int step = 0; step < matrixSize - 1; ++step)
    {
        // Calculate w(step).
        double normal = FindNormal (A, matrixSize, step);
        for (int index = 0; index < matrixSize; ++index)
        {
            if (index < step)
            {
                w[index] = 0.0;
            }
            else
            {
                w[index] = index == step ? A[step][step] - normal : A[index][step];
            }
        }

        normal = FindVectorNormal (w, matrixSize);
        for (int index = 0; index < matrixSize; ++index)
        {
            w[index] /= normal;
        }

        // Apply w(step) to matrices.
        for (int col = step; col < matrixSize; ++col)
        {
            ApplyW (A, matrixSize, col, step, w);
        }

        for (int col = 0; col < results; ++col)
        {
            ApplyW (B, matrixSize, col, step, w);
        }
    }

    for (int step = matrixSize - 1; step >= 0; --step)
    {
        MultiplyRow (B, matrixSize, results, step, 1.0 / A[step][step]);
        for (int row = step - 1; row >= 0; --row)
        {
            AddMultipliedRow (B, matrixSize, results, row, step, -A[row][step]);
        }
    }

    return true;
}
