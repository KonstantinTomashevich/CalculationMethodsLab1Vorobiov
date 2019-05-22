#include "matrixutils.h"
#include "mtwister.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern MTRand *GlobalRand;
static double MinValue = 0.0;
static double MaxValue = 0.0;
static double Addition = 0.0;
static double Numberer = 0.0;

double RandomMatrixValue ()
{
    if (MinValue == MaxValue)
    {
        MaxValue = pow (2.0, N / 4.0);
        MinValue = -MaxValue;
        Addition = pow (10.0, -13.0);
        Numberer = pow (10.0, 13.0);
    }

    double result = GenRand (GlobalRand) * (MaxValue - MinValue) + MinValue;
    double asNumber;

    if (modf (result * Numberer, &asNumber) < 0.0001)
    {
        result += Addition;
    }

    return result;
}

double **AllocateMatrix (int rows, int cols)
{
    double **matrix = calloc (rows, sizeof (double *));
    for (int index = 0; index < rows; ++index)
    {
        matrix[index] = calloc (cols, sizeof (double));
    }

    return matrix;
}

double **CopyMatrix (double **matrix, int rows, int cols)
{
    double **copy = AllocateMatrix (rows, cols);
    CopyMatrixInto (matrix, rows, cols, copy);
    return copy;
}

void CopyMatrixInto (double **matrix, int rows, int cols, double **output)
{
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            output[row][col] = matrix[row][col];
        }
    }
}

double **TransformMatrixByRowOrder (double **matrix, int rows, int cols, int *rowOrder)
{
    double **newMatrix = AllocateMatrix (rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            newMatrix[row][col] = matrix[rowOrder[row]][col];
        }
    }

    return newMatrix;
}

double **TransformMatrixByColOrder (double **matrix, int rows, int cols, int *colOrder)
{
    double **newMatrix = AllocateMatrix (rows, cols);
    for (int col = 0; col < cols; ++col)
    {
        for (int row = 0; row < rows; ++row)
        {
            newMatrix[row][col] = matrix[row][colOrder[col]];
        }
    }

    return newMatrix;
}

double **TransposeMatrix (double **matrix, int rows, int cols)
{
    double **transposed = AllocateMatrix (cols, rows);
    for (int row = 0; row < cols; ++row)
    {
        for (int col = 0; col < rows; ++col)
        {
            transposed[row][col] = matrix[col][row];
        }
    }

    return transposed;
}

void FreeMatrix (double **matrix, int rows, int cols)
{
    for (int row = 0; row < rows; ++row)
    {
        free (matrix[row]);
    }

    free (matrix);
}

double ColumnQuadricNormal (double **matrix, int rows, int colIndex)
{
    double result = 0.0;
    for (int row = 0; row < rows; ++row)
    {
        result += matrix[row][colIndex] * matrix[row][colIndex];
    }

    return sqrt (result);
}

double ColumnsScalarMultiplication (double **matrix, int rows, int col1Index, int col2Index)
{
    double result = 0.0;
    for (int row = 0; row < rows; ++row)
    {
        result += matrix[row][col1Index] * matrix[row][col2Index];
    }

    return result;
}

void FillTaskSpecificMatrix (double **matrix)
{
    for (int col = 1; col < MATRIX_SIZE; ++col)
    {
        for (int row = 0; row < col; ++row)
        {
            double value = RandomMatrixValue ();
            matrix[row][col] = value;
            matrix[col][row] = value;
        }
    }

    for (int diag = 0; diag < MATRIX_SIZE; ++diag)
    {
        double value = 0.0;
        for (int col = 0; col < MATRIX_SIZE; ++col)
        {
            if (col != diag)
            {
                value += fabs (matrix[diag][col]);
            }
        }

        matrix[diag][diag] = value;
    }
}

void FillDefaultMatrix (double **matrix, int rows, int cols)
{
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            matrix[row][col] = RandomMatrixValue ();
        }
    }
}

void PrintMatrix (double **matrix, int rows, int cols)
{
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            printf ("%20.13lf ", matrix[row][col]);
        }

        printf ("\n");
    }
}

void MultiplyRow (double **matrix, int rows, int cols, int row, double by)
{
    MultiplyRowPart (matrix, rows, cols, row, by, 0);
}

void MultiplyRowPart (double **matrix, int rows, int cols, int row, double by, int start)
{
    for (int col = start; col < cols; ++col)
    {
        matrix[row][col] *= by;
    }
}

void MultiplyCol (double **matrix, int rows, int cols, int col, double by)
{
    for (int row = 0; row < rows; ++row)
    {
        matrix[row][col] *= by;
    }
}

void SwapRows (double **matrix, int rows, int cols, int row1, int row2)
{
    double *temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
}

void SwapCols (double **matrix, int rows, int cols, int col1, int col2)
{
    for (int row = 0; row < rows; ++row)
    {
        double temp = matrix[row][col1];
        matrix[row][col1] = matrix[row][col2];
        matrix[row][col2] = temp;
    }
}

void MultiplyMatrices (double **first, double **second, double **output, int firstRows, int firstCols, int secondsCols)
{
    for (int row = 0; row < firstRows; ++row)
    {
        for (int col = 0; col < secondsCols; ++col)
        {
            double value = 0.0;
            for (int iterator = 0; iterator < firstCols; ++iterator)
            {
                value += first[row][iterator] * second[iterator][col];
            }

            output[row][col] = value;
        }
    }
}

void AddMultipliedRow (double **matrix, int rows, int cols, int dst, int src, double modifier)
{
    for (int col = 0; col < cols; ++col)
    {
        matrix[dst][col] += matrix[src][col] * modifier;
    }
}

void AddMultipliedRowPart (double **matrix, int rows, int cols, int dst, int src, double modifier, int startFrom)
{
    for (int col = startFrom; col < cols; ++col)
    {
        matrix[dst][col] += matrix[src][col] * modifier;
    }
}

void AddMultipliedCol (double **matrix, int rows, int cols, int dst, int src, double modifier)
{
    for (int row = 0; row < rows; ++row)
    {
        matrix[row][dst] += matrix[row][src] * modifier;
    }
}
