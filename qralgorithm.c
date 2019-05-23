#include "qralgorithm.h"
#include "matrixutils.h"
#include "def.h"

#include <math.h>
#include <stdlib.h>

#define BARRIER (1.0/10000)
#define DEFLATION_BARRIER (1.0/1000000000000)
#define INCORRECT_ROTATION -2

void SimilarTransformation (double **A, int row, int col, int matrixSize)
{
    double modifier = A[row][col] / A[row][row - 1];
    AddMultipliedCol (A, matrixSize, matrixSize, col, row - 1, -modifier);
    AddMultipliedRow (A, matrixSize, matrixSize, row - 1, col, modifier);
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
            if (row == index)
            {
                value += y * copy[index][col];
                value += x * copy[index + 1][col];
            }
            else if (row == index + 1)
            {
                value -= x * copy[index][col];
                value += y * copy[index + 1][col];
            }
            else
            {
                value = copy[row][col];
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
            if (col == index)
            {
                value += copy[row][index] * y;
                value += copy[row][index + 1] * x;
            }
            else if (col == index + 1)
            {
                value -= copy[row][index] * x;
                value += copy[row][index + 1] * y;
            }
            else
            {
                value = copy[row][col];
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
    double *rotationsBuffer = calloc ((matrixSize - 1) * 2, sizeof (double));
    double *previousDiag = calloc (matrixSize, sizeof (double));

    while (iteration < pow (matrixSize, 3))
    {
        for (int index = 0; index < matrixSize - 1; ++index)
        {
            RotationTransformation (A, matrixSize, index, rotationsBuffer);
        }

        for (int index = 0; index < matrixSize - 1; ++index)
        {
            UndoRotationTransformation (A, matrixSize, index, rotationsBuffer);
        }

        int found = 0;
        int complex = 0;

        for (int testIndex = 0; testIndex < matrixSize; ++testIndex)
        {
            if (testIndex + 1 < matrixSize && fabs(A[testIndex + 1][testIndex]) > DEFLATION_BARRIER)
            {
                found += 2;
                complex += 2;
                ++testIndex;
                continue;
            }

            for (int searchIndex = 0; searchIndex < matrixSize; ++searchIndex)
            {
                if (fabs (A[testIndex][testIndex] - previousDiag[searchIndex]) < BARRIER)
                {
                    found++;
                    break;
                }
            }
        }

        if (found == matrixSize && found - complex > 0)
        {
            break;
        }

        for (int index = 0; index < matrixSize; ++index)
        {
            previousDiag[index] = A[index][index];
        }

        ++iteration;
    }

    free (previousDiag);
}
