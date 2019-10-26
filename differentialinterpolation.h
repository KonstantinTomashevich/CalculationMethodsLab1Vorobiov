#pragma once

void DifferentialInterpolation (double (*F) (double), double min, double max, double *X, int points, double **coeffs);
void ChebyshevDifferentialInterpolation (double (*F) (double), double min, double max, int points, double **coeffs);
