#pragma once
// Output is [power+1]x1 matrix (because Householder is used).
void QuadraticRegression(const double *x, const double *y, int count, int power, double ***output);
