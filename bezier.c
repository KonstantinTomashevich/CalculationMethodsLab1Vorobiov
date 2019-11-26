#include "bezier.h"
#include <stdlib.h>

/// output (t) = (1 - t) * first (t) + t * second (t)
static void Sum (const double *first, const double *second, int size, double *output)
{
    output[0] = first[0];
    for (int index = 1; index < size; ++index)
    {
        output[index] = first[index] + second[index - 1] - first[index - 1];
    }
}

static void BuildForAxis (const double *values, int count, double *output)
{
    double *matrix = calloc (count * count, sizeof (double));

    for (int index = 0; index < count; ++index)
    {
        matrix[index * count] = values[index];
    }

    for (int rows = count; rows > 1; --rows)
    {
        for (int row = 1; row < rows; ++row)
        {
            Sum (matrix + (row - 1) * count, matrix + row * count, count, output);
            if (rows > 2)
            {
                for (int col = 0; col < count; ++col)
                {
                    matrix[(row - 1) * count + col] = output[col];
                }
            }
        }
    }

    free (matrix);
}

void BezierCasteljau (const double *x, const double *y, int count, double **xOut, double **yOut)
{
    *xOut = calloc (count, sizeof (double));
    *yOut = calloc (count, sizeof (double));
    BuildForAxis (x, count, *xOut);
    BuildForAxis (y, count, *yOut);
}
