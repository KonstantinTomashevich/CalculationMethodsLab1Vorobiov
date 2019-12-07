#pragma once

double RungeRule (double (*Integrator) (double (*) (double), int, double, double, long),
    double (*F) (double), int power, double a, double b, double epsilon, unsigned long *resultParts);
