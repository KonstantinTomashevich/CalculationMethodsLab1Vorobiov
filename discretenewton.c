#include "discretenewton.h"
#include <math.h>

#define H (1.0/10000)
#define BARRIER (3.0/10000000000000000)

bool DiscreteNewton (double (*F) (double), double *x, double a, double b, int *iterations)
{
    *iterations = 1;
    *x = b;

    while (true)
    {
        double newX = *x - H * F (*x) / (F (*x + H) - F (*x));
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
