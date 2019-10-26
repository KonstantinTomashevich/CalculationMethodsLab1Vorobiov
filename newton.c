#include "discretenewton.h"
#include <math.h>

#define BARRIER (3.0/10000000000000000)

bool Newton (double (*F) (double), double (*dF) (double), double *x, double a, double b, int *iterations)
{
    *iterations = 1;
    while (true)
    {
        double newX = *x - F (*x) / dF (*x);
        if (newX < a || newX > b)
        {
            return false;
        }

        if (fabs (newX - *x) < BARRIER)
        {
            *x = newX;
            return true;
        }

        *x = newX;
        ++(*iterations);
    }
}
