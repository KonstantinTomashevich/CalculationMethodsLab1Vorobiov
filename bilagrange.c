#include "bilagrange.h"

static void PrintLambdaAndCalculateMultiplier (double *values, int count, int index, const char *var, FILE *output)
{
    for (int i = 0; i < count; ++i)
    {
        if (index != i)
        {
            fprintf (output, "*((%s-%0.16lf)/%lf)", var, values[i], values[index] - values[i]);
        }
    }
}

void PrintBiLagrangePolynomial (double *x, double *y, int gridSize, double **z, FILE *output)
{
    for (int row = 0; row < gridSize; ++row)
    {
        for (int col = 0; col < gridSize; ++col)
        {
            fprintf (output, "%0.16lf", z[row][col]);
            PrintLambdaAndCalculateMultiplier (x, gridSize, row, "x", output);
            PrintLambdaAndCalculateMultiplier (y, gridSize, col, "y", output);

            if (row < gridSize - 1 || col < gridSize - 1)
            {
                fprintf (output, " + ");
            }
        }
    }
}
