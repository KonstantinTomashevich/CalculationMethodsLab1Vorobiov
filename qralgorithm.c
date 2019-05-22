#include "qralgorithm.h"
#include "matrixutils.h"
#include "def.h"

#include <math.h>
#include <stdlib.h>

#define BARRIER 0.001
#define DEFLATION_BARRIER (1.0/10000000000000)
#define INCORRECT_ROTATION -2

void SimilarTransformation (double **A, int row, int col, int matrixSize)
{
    AddMultipliedCol (A, matrixSize, matrixSize, col, row - 1, -A[row][col] / A[row][row - 1]);
    AddMultipliedRow (A, matrixSize, matrixSize, row - 1, col, A[row][col] / A[row][row - 1]);
}

void RotationTransformation (double **A, int matrixSize, int index, double *undoBuffer)
{
    double a = A[index][index];
    double b = A[index + 1][index];

    if (fabs (b) < DEFLATION_BARRIER)
    {
        undoBuffer[index * 2] = INCORRECT_ROTATION;
        undoBuffer[index * 2 + 1] = INCORRECT_ROTATION;
        return;
    }

    double y = sqrt (a * a + b * b) / (a + b * b / a);
    double x = b * y / a;

    undoBuffer[index * 2] = x;
    undoBuffer[index * 2 + 1] = y;

    double **copy = CopyMatrix (A, matrixSize, matrixSize);
    for (int row = 0; row < matrixSize; ++row)
    {
        for (int col = 0; col < matrixSize; ++col)
        {
            double value = 0.0;
            for (int iterator = 0; iterator < matrixSize; ++iterator)
            {
                double firstValue;
                if ((row == index && iterator == index) ||
                    (row == index + 1 && iterator == index + 1))
                { firstValue = y; }
                else if (row == index && iterator == index + 1)
                { firstValue = x; }
                else if (row == index + 1 && iterator == index)
                { firstValue = -x; }
                else if (row == iterator)
                { firstValue = 1.0; }
                else
                { firstValue = 0.0; }

                value += firstValue * copy[iterator][col];
            }

            A[row][col] = value;
        }
    }

    FreeMatrix (copy, matrixSize, matrixSize);
}

void UndoRotationTransformation (double **A, int matrixSize, int index, double *undoBuffer)
{
    double x = undoBuffer[index * 2];
    double y = undoBuffer[index * 2 + 1];

    if (x == INCORRECT_ROTATION || y == INCORRECT_ROTATION)
    {
        return;
    }
    
    double **copy = CopyMatrix (A, matrixSize, matrixSize);
    for (int row = 0; row < matrixSize; ++row)
    {
        for (int col = 0; col < matrixSize; ++col)
        {
            double value = 0.0;
            for (int iterator = 0; iterator < matrixSize; ++iterator)
            {
                double secondValue;
                if ((iterator == index && col == index) ||
                    (iterator == index + 1 && col == index + 1))
                { secondValue = y; }
                else if (iterator == index && col == index + 1)
                { secondValue = -x; }
                else if (iterator == index + 1 && col == index)
                { secondValue = x; }
                else if (iterator == col)
                { secondValue = 1.0; }
                else
                { secondValue = 0.0; }

                value += copy[row][iterator] * secondValue;
            }

            A[row][col] = value;
        }
    }

    FreeMatrix (copy, matrixSize, matrixSize);
}

void QRAlgorithm (double **A, int matrixSize)
{
    for (int row = matrixSize - 1; row >= 2; --row)
    {
        for (int col = row - 2; col >= 0; --col)
        {
            SimilarTransformation (A, row, col, matrixSize);
        }
    }

    int iteration = 0;
    double previousNormal = 0.0;
    double *rotationsBuffer = calloc ((matrixSize - 1) * 2, sizeof (double));

    while (iteration < pow (matrixSize, 4))
    {
        for (int index = 0; index < matrixSize - 1; ++index)
        {
            RotationTransformation (A, matrixSize, index, rotationsBuffer);
        }

        for (int index = 0; index < matrixSize - 1; ++index)
        {
            UndoRotationTransformation (A, matrixSize, index, rotationsBuffer);
        }

        double normal = 0.0;
        for (int index = 0; index < matrixSize; ++index)
        {
            normal += pow (A[index][index], 2);
        }

        normal = sqrt (normal);
        if (iteration > 0 && fabs (normal - previousNormal) < BARRIER)
        {
            break;
        }

        previousNormal = normal;
        ++iteration;
    }
}
