#pragma once

void CubicSpline (double (*F) (double), double min, double max, int splineCount, double **output);
void CubicSplineForEquidistantPoints (const double *y, double step, int splineCount, double **output);
