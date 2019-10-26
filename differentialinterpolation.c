#include "differentialinterpolation.h"
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void DifferentialInterpolation (double (*F) (double), double min, double max, int points, double **coeffs)
{
    double *x = malloc (sizeof (double) * points);
    double *y = malloc (sizeof (double) * points * points);
    *coeffs = malloc (sizeof (double) * points);

    double step = (max - min) / (points - 1);
    for (int index = 0; index < points; ++index)
    {
        x[index] = min + step * index;
        y[index] = F (x[index]);
    }

    for (int col = 1; col < points; ++col)
    {
        for (int row = 0; row < points - col; ++row)
        {
            int i1 = row;
            int i2 = i1 + col;
            int sourceCol = col - 1;
            y[col * points + row] =
                (y[sourceCol * points + i2 - sourceCol] - y[sourceCol * points + i1]) / (x[i2] -x[i1]);
        }
    }

    for (int index = 0; index < points; ++index)
    {
        (*coeffs)[index] = y[index * points];
    }

    free (x);
    free (y);
}
