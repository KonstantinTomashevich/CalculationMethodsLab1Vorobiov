#pragma once

void BiCubicSpline (double (*G) (double, double), double minX, double maxX, double minY, double maxY,
    int gridSize, double ***xOutput);
double BiCubicSplineCalculate (double x, double y, double minX, double maxX, double minY, double maxY,
    int gridSize, double **xResult);
