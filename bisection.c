#include "bisection.h"
#include <math.h>
#define BARRIER 0.0001

static int sign (double x)
{
    return x < 0 ? -1 : (x == 0 ? 0 : 1);
}

void Bisection (double (*F) (double), double *a, double *b, int *iteration)
{
    double left = F (*a);
    double right = F (*b);
    *iteration = 0;

    while (fabs (*a - *b) > BARRIER)
    {
        double middle = *a + fabs (*a - *b) / 2;
        double value = F (middle);
        int valueSign = sign (value);

        if (valueSign == 0)
        {
            break;
        }
        else if (valueSign == sign (left))
        {
            left = value;
            *a = middle;
        }
        else if (valueSign == sign (right))
        {
            right = value;
            *b = middle;
        }

        ++*iteration;
    }
}
