#pragma once
#include <stdbool.h>

bool Newton (double (*F) (double), double (*dF) (double), double *x, double a, double b, int *iterations);
