#include "runge.h"
#include <math.h>

double RungeRule (double (*Integrator) (double (*) (double), int, double, double, long),
                  double (*F) (double), int power, double a, double b, double epsilon, unsigned long *resultParts)
{
    double rMultiplier = 1.0 / (pow (2.0, power * 2 - 1) - 1);
    unsigned long parts = 1;
    double previous;
    double next = Integrator (F, power, a, b, parts);

    do
    {
        previous = next;
        parts = parts << 1u;
        next = Integrator (F, power, a, b, parts);
    } while (fabs (next - previous) * rMultiplier > epsilon);

    *resultParts = parts;
    return next;
}
