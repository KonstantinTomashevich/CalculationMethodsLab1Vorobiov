#include "runge.h"
#include <math.h>
#include <stdio.h>

double RungeRule (double (*Integrator) (double (*) (double), int, double, double, long),
                  double (*F) (double), int power, int accuracy, double a, double b, double epsilon,
                  unsigned long *resultParts, unsigned long *rungeSteps)
{
    char firstRun = 1;
    unsigned long firstParts = 1;
    unsigned long previousFirstParts = 1;
    unsigned long steps = 1;

    double first;
    double second;

    double previousC;
    double c;
    double h;

    do
    {
        if (!firstRun)
        {
            ++steps;
            double nextH = pow (epsilon / c, 1.0 / accuracy);

            if (!isfinite (nextH))
            {
                printf ("Error, next step became too small! Returning last correct result...\n");
                *rungeSteps = steps - 1;
                return first;
            }

            previousFirstParts = firstParts;
            firstParts = lround ((b - a) / nextH);
        }

        first = Integrator (F, power, a, b, firstParts);
        second = Integrator (F, power, a, b, firstParts * 2);
        h = (b - a) / firstParts;

        previousC = c;
        c = fabs ((first - second) / (pow (h / 2, accuracy) - pow (h, accuracy)));

        if (!firstRun && c * pow (h, accuracy) > previousC * pow ((b - a) / previousFirstParts, accuracy))
        {
            printf ("Found accuracy degradation! Returning last correct result...\n");
            *resultParts = previousFirstParts;
            *rungeSteps = steps;
            return Integrator (F, power, a, b, previousFirstParts);
        }

        firstRun = 0;
    } while (c * pow (h, accuracy) > epsilon);

    *resultParts = firstParts;
    *rungeSteps = steps;
    return first;
}
