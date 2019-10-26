#include "differentialinterpolation.h"
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void DifferentialInterpolation (double (*F) (double), double min, double max, int points, double **coeffs)
{
    double *x = malloc (sizeof (double) * points);
    double *y = malloc (sizeof (double) * (points + 1) * points / 2);
    *coeffs = malloc (sizeof (double) * points);

    double step = (max - min) / (points - 1);
    for (int index = 0; index < points; ++index)
    {
        x[index] = min + step * index;
        y[index] = F (x[index]);
    }

    for (int col = 1; col < points; ++col)
    {
        int sourceCol = col - 1;
        int sourceColOffset = (points * 2 - sourceCol + 1) * sourceCol / 2;
        int targetColOffset = (points * 2 - col + 1) * col / 2;

        for (int row = 0; row < points - col; ++row)
        {
            int i1 = row;
            int i2 = i1 + col;
            y[targetColOffset + row] = (y[sourceColOffset + i2 - sourceCol] - y[sourceColOffset + i1]) / (x[i2] - x[i1]);
        }
    }


    for (int index = 0, tableIndex = 0; index < points; tableIndex += points - index, ++index)
    {
        (*coeffs)[index] = y[tableIndex];
    }

    free (x);
    free (y);
}
