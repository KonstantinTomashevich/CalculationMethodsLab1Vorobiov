#include "bicubicspline.h"
#include "matrixutils.h"
#include "cubicspline.h"
#include <stdlib.h>
#include <math.h>

void BiCubicSpline (double (*G) (double, double), double minX, double maxX, double minY, double maxY,
                    int gridSize, double ***xOutput)
{
    double **zMatrix = AllocateMatrix (gridSize, gridSize);
    double x = minX;
    double y = minY;

    double stepX = (maxX - minX) / (gridSize - 1);
    double stepY = (maxY - minY) / (gridSize - 1);

    for (int row = 0; row < gridSize; ++row)
    {
        for (int col = 0; col < gridSize; ++col)
        {
            zMatrix[row][col] = G (x, y);
            y += stepY;
        }

        y = minY;
        x += stepX;
    }

    double **ys = calloc (gridSize, sizeof (double *));
    for (int row = 0; row < gridSize; ++row)
    {
        CubicSplineForEquidistantPoints (zMatrix[row], stepY, gridSize - 1, ys + row);
    }

    double *nextSplinesValues = calloc (gridSize, sizeof (double));
    *xOutput = calloc ((gridSize - 1) * 4, sizeof (double *));

    for (int row = 0; row < (gridSize - 1) * 4; ++row)
    {
        for (int paramIndex = 0; paramIndex < gridSize; ++paramIndex)
        {
            nextSplinesValues[paramIndex] = ys[paramIndex][row];
        }

        CubicSplineForEquidistantPoints (nextSplinesValues, stepX, gridSize - 1, (*xOutput) + row);
    }

    FreeMatrix (zMatrix, gridSize, gridSize);
    free (nextSplinesValues);

    for (int row = 0; row < gridSize; ++row)
    {
        free (ys[row]);
    }

    free (ys);
}

static double SplineResult (double *coefs, double value, double pivot)
{
    return coefs[0] + coefs[1] * (value - pivot) +
        coefs[2] * pow (value - pivot, 2) +
        coefs[3] * pow (value - pivot, 3);
}

double BiCubicSplineCalculate (double x, double y, double minX, double maxX, double minY, double maxY,
                               int gridSize, double **xResult)
{
    if (x < minX || x > maxX || y < minY || y > maxY)
    {
        return 0.0;
    }

    int xSegmentIndex = (int) trunc(((x - minX) / (maxX - minX)) * (gridSize - 1));
    if (xSegmentIndex == gridSize - 1)
    {
        xSegmentIndex = gridSize - 2;
    }

    int ySegmentIndex = (int) trunc(((y - minY) / (maxY - minY)) * (gridSize - 1));
    if (ySegmentIndex == gridSize - 1)
    {
        ySegmentIndex = gridSize - 2;
    }

    double pivotX = minX + (maxX - minX) * xSegmentIndex / (gridSize - 1);
    double a = SplineResult (xResult[ySegmentIndex * 4] + xSegmentIndex * 4, x, pivotX);
    double b = SplineResult (xResult[ySegmentIndex * 4 + 1] + xSegmentIndex * 4, x, pivotX);
    double c = SplineResult (xResult[ySegmentIndex * 4 + 2] + xSegmentIndex * 4, x, pivotX);
    double d = SplineResult (xResult[ySegmentIndex * 4 + 3] + xSegmentIndex * 4, x, pivotX);

    double pivotY = minY + (maxY - minY) * ySegmentIndex / (gridSize - 1);
    return a + b * (y - pivotY) + c * pow (y - pivotY, 2) + d * pow (y - pivotY, 3);
}
