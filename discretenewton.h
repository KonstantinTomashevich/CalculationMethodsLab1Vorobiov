#pragma once
#include <stdbool.h>

bool DescriteNewton (double (*F) (double), double *x, double a, double b, int *iterations);
