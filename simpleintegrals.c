#include "simpleintegrals.h"
#include <stdio.h>

double LeftRectangleIntegral (double (*F) (double), double a, double b, long parts)
{
    double x = a;
    double sum = 0.0;
    double step = (b - a) / parts;

    while (parts--)
    {
        sum += F (x) * step;
        x += step;
    }

    return sum;
}

double RightRectangleIntegral (double (*F) (double), double a, double b, long parts)
{
    double step = (b - a) / parts;
    double x = a + step;
    double sum = 0.0;

    while (parts--)
    {
        sum += F (x) * step;
        x += step;
    }

    return sum;
}

double MediumRectangleIntegral (double (*F) (double), double a, double b, long parts)
{
    double step = (b - a) / parts;
    double x = a + step / 2.0;
    double sum = 0.0;

    while (parts--)
    {
        sum += F (x) * step;
        x += step;
    }

    return sum;
}

double TrapeziumRectangleIntegral (double (*F) (double), double a, double b, long parts)
{
    double x = a;
    double sum = 0.0;
    double step = (b - a) / parts;

    while (parts--)
    {
        sum += (F (x) + F (x + step)) * 0.5 * step;
        x += step;
    }

    return sum;
}

double SimpsonRectangleIntegral (double (*F) (double), double a, double b, long parts)
{
    double x = a;
    double sum = 0.0;
    double step = (b - a) / parts;

    while (parts--)
    {
        sum += (F (x) + F (x + step) + 4.0 * F (x + step / 2.0)) * step / 6.0;
        x += step;
    }

    return sum;
}

static double newton5[] = {0.5 * 7 / 45, 0.5 * 32 / 45, 0.5 * 12 / 45, 0.5 * 32 / 45, 0.5 * 7 / 45};

static double newton7[] =
    {0.5 * 41 / 420, 0.5 * 216 / 420, 0.5 * 27 / 420, 0.5 * 272 / 420, 0.5 * 27 / 420, 0.5 * 216 / 420, 0.5 * 41 / 420};

static double newton9[] = {0.5 * 989 / 14175,
                           0.5 * 5888 / 14175,
                           0.5 * -928 / 14175,
                           0.5 * 10496 / 14175,
                           0.5 * -4540 / 14175,
                           0.5 * 10496 / 14175,
                           0.5 * -928 / 14175,
                           0.5 * 5888 / 14175,
                           0.5 * 989 / 14175};

static double newton11[] = {0.5 * 0.05366829673142116,
                            0.5 * 0.3550718828218748,
                            0.5 * -0.1620871411834986,
                            0.5 * 0.9098925764163387,
                            0.5 * -0.8703102450914472,
                            0.5 * 1.4275292606105543,
                            0.5 * -0.8703102450914436,
                            0.5 * 0.9098925764163849,
                            0.5 * -0.1620871411834854,
                            0.5 * 0.3550718828218774,
                            0.5 * 0.05366829673142342};

static double newton13[] = {0.5 * 0.043278974999896684,
                            0.5 * 0.31407221347352104,
                            0.5 * -0.24064392741857343,
                            0.5 * 1.1329977956512114,
                            0.5 * -1.633011274053057,
                            0.5 * 2.7755193372391824,
                            0.5 * -2.7844262397844264,
                            0.5 * 2.7755193372392784,
                            0.5 * -1.6330112740530711,
                            0.5 * 1.1329977956512016,
                            0.5 * -0.24064392741856697,
                            0.5 * 0.3140722134735046,
                            0.5 * 0.043278974999898155};

static double newton15[] = {0.5 * 0.0360689424375301,
                            0.5 * 0.28417558935857096,
                            0.5 * -0.30805069400994944,
                            0.5 * 1.399497820605065,
                            0.5 * -2.647995210668152,
                            0.5 * 5.048155507760063,
                            0.5 * -6.7157289774485776,
                            0.5 * 7.8077540439305,
                            0.5 * -6.715728977448385,
                            0.5 * 5.0481555077600735,
                            0.5 * -2.6479952106679328,
                            0.5 * 1.3994978206050246,
                            0.5 * -0.30805069400996593,
                            0.5 * 0.2841755893585973,
                            0.5 * 0.0360689424375374};

static double *newtonCoefficients[] = {newton5, newton7, newton9, newton11, newton13, newton15};

double NewtonIntegral (double (*F) (double), int power, double a, double b, long parts)
{
    if (power < 5 || power % 2 == 0 || power > 15)
    {
        printf ("Newton integral power must be >= 5, <= 15 and odd!\n");
        return 0.0;
    }

    double x = a;
    double sum = 0.0;
    double step = (b - a) / parts;
    double partialStep = step / (power - 1);
    double *coefficients = newtonCoefficients[(power - 5) / 2];

    while (parts--)
    {
        for (int index = 0; index < power; ++index)
        {
            sum += F (x) * coefficients[index] * step;
            x += partialStep;
        }

        x -= partialStep;
    }

    return sum;
}
