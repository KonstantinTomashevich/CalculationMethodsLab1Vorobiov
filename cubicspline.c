#include "cubicspline.h"
#include <stdlib.h>
#include <math.h>

void CubicSpline (double (*F) (double), double min, double max, int splineCount, double **output)
{
    double *x = malloc (sizeof (double) * (splineCount + 1));
    double *y = malloc (sizeof (double) * (splineCount + 1));
    *output = malloc (sizeof (double) * splineCount * 4); // for each spline -- (a, b, c, d)

    double step = (max - min) / splineCount;
    for (int index = 0; index <= splineCount; ++index)
    {
        x[index] = min + step * index;
        y[index] = F (x[index]);
    }

    double *d = malloc (sizeof (double) * (splineCount + 1));
    double *b = malloc (sizeof (double) * (splineCount + 1));

    d[0] = 1;
    b[0] = 0;

    for (int index = 1; index < splineCount; ++index)
    {
        d[index] = 2;
        b[index] = 1.5 * ((y[index + 1] - y[index]) / step - (y[index] - y[index - 1]) / step) / step;
    }

    d[splineCount] = 1;
    b[splineCount] = 0;

    for (int index = 2; index < splineCount; ++index)
    {
        double multiplier = 0.5 / d[index - 1];
        d[index] -= 0.5 * multiplier;
        b[index] -= b[index - 1] * multiplier;
    }

    (*output)[4 * (splineCount - 1) + 2] = 0.0;
    (*output)[4 * (splineCount - 2) + 2] = b[splineCount - 1] / d[splineCount - 1];

    for (int index = splineCount - 2; index >= 1; --index)
    {
        (*output)[4 * (index - 1) + 2] = (b[index] - 0.5 * (*output)[4 * index + 2]) / d[index];
    }

    for (int index = 0; index < splineCount; ++index)
    {
        (*output)[4 * index] = y[index];
        if (index == 0)
        {
            (*output)[4 * index + 3] = (*output)[4 * index + 2] / step / 3.0;
            (*output)[4 * index + 1] = (y[index + 1] - y[index]) / step -
                step * (*output)[4 * index + 2] / 3.0;
        }
        else
        {
            (*output)[4 * index + 3] = ((*output)[4 * index + 2] - (*output)[4 * (index - 1) + 2]) / step / 3.0;
            (*output)[4 * index + 1] = (y[index + 1] - y[index]) / step -
                step * ((*output)[4 * index + 2] + 2.0 * (*output)[4 * (index - 1) + 2]) / 3.0;
        }
    }

    for (int index = splineCount - 1; index > 0; --index)
    {
        (*output)[4 * (index) + 2] = (*output)[4 * (index - 1) + 2];
    }

    (*output)[2] = 0.0;

    free (x);
    free (y);
    free (d);
    free (b);
}
