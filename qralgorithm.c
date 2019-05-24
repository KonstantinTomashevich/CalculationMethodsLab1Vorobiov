#include "qralgorithm.h"
#include "matrixutils.h"
#include "def.h"

#include <math.h>
#include <stdlib.h>

#define BARRIER (1.0/10000)
#define DEFLATION_BARRIER (1.0/1000000000)
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

int QRAlgorithm (double **A, int matrixSize)
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
    double *previousValues = calloc (matrixSize * 2, sizeof (double));
    double *currentValues = calloc (matrixSize * 2, sizeof (double));

    while (1)
    {
        for (int index = 0; index < matrixSize - 1; ++index)
        {
            RotationTransformation (A, matrixSize, index, rotationsBuffer);
        }

        for (int index = 0; index < matrixSize - 1; ++index)
        {
            UndoRotationTransformation (A, matrixSize, index, rotationsBuffer);
        }

        int checkIndex = 0;
        while (checkIndex < matrixSize)
        {
            double realValue;
            double complexValue;

            if (checkIndex < matrixSize - 1 && fabs (A[checkIndex + 1][checkIndex]) > DEFLATION_BARRIER)
            {
                realValue = (A[checkIndex][checkIndex] + A[checkIndex + 1][checkIndex + 1]) / 2;
                complexValue = sqrt (fabs (
                    pow (A[checkIndex][checkIndex] - A[checkIndex + 1][checkIndex + 1], 2) +
                        4 * A[checkIndex][checkIndex + 1] * A[checkIndex + 1][checkIndex])) / 2;

                currentValues[checkIndex * 2] = realValue;
                currentValues[checkIndex * 2 + 1] = complexValue;
                ++checkIndex;

                currentValues[checkIndex * 2] = realValue;
                currentValues[checkIndex * 2 + 1] = -complexValue;
                ++checkIndex;
            }
            else
            {
                realValue = A[checkIndex][checkIndex];
                complexValue = 0.0;

                currentValues[checkIndex * 2] = realValue;
                currentValues[checkIndex * 2 + 1] = complexValue;
                ++checkIndex;
            }
        }

        int found = 0;
        for (int testIndex = 0; testIndex < matrixSize; ++testIndex)
        {
            double real = currentValues[testIndex * 2];
            double complex = currentValues[testIndex * 2 + 1];

            for (int searchIndex = 0; searchIndex < matrixSize; ++searchIndex)
            {
                if (fabs (real - previousValues[searchIndex * 2]) < BARRIER &&
                    fabs (complex - previousValues[searchIndex * 2 + 1]) < BARRIER)
                {
                    ++found;
                    break;
                }
            }

            if (found <= testIndex)
            {
                break;
            }
        }

        if (found == matrixSize)
        {
            break;
        }

        for (int index = 0; index < matrixSize * 2; ++index)
        {
            previousValues[index] = currentValues[index];
        }

        ++iteration;
    }

    free (previousValues);
    free (currentValues);
    return iteration;
}
