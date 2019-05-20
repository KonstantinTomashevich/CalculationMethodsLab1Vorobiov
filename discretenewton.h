#pragma once
#include <stdbool.h>

bool DiscreteNewton (double (*F) (double), double *x, double a, double b, int *iterations);
