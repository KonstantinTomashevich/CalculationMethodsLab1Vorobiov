#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "bisection.h"
#include "discretenewton.h"
#include "newton.h"
#include "differentialinterpolation.h"
#include "cubicspline.h"

#define SEGMENT_COUNT 3

#define INTERPOLATION_NODE_SET_COUNT 5

#define INTERPOLATION_RESULT_MUTE_AFTER 2

static const char *pythonSupportCode =
    "import numpy as np\n"
    "import math\n"
    "import random as rnd\n"
    "from matplotlib import pyplot as plt\n"
    "import scipy.linalg as lg\n"
    "import time\n"
    "import copy\n"
    "\n"
    "def f(x):\n"
    "    return (x**9 + math.pi)*math.cos(math.log((x**2)+1))/math.exp(x**2) - x/2019\n"
    "\n"
    "def show_plot(Xs, Ys, xlims=[-4,4], ylims=[-4,4], show_f = False):\n"
    "    plt.plot(Xs, Ys)\n"
    "    plt.axvline(0, color='black')\n"
    "    plt.axhline(0, color='black')\n"
    "    plt.grid()\n"
    "    plt.ylim(ylims)\n"
    "    plt.xlim(xlims)\n"
    "    if show_f:\n"
    "        plt.plot(Xs, [f(x) for x in Xs])\n"
    "    plt.show()\n"
    "\n"
    "def find_spline_index(x, t):\n"
    "    n = len(x)\n"
    "    k = 0\n"
    "    while k < n:\n"
    "        mid = (k + n) // 2\n"
    "        if t < x[mid]:\n"
    "            n = mid\n"
    "        else:\n"
    "            k = mid + 1\n"
    "    return k\n"
    "\n"
    "def build_spline(a, b, c, d):\n"
    "    def spline(x,t):\n"
    "        i = find_spline_index(x, t) - 1\n"
    "        dx = t - x[i]\n"
    "        result = a[i] + b[i] * dx + c[i] * dx ** 2.0 + d[i] * dx ** 3.0\n"
    "        return result\n"
    "    return spline\n"
    "\n"
    "def generate_xs(left, right, n):\n"
    "    return sorted([rnd.uniform(left,right) for i in range(n)])\n";

int interpolationNodeSet[INTERPOLATION_NODE_SET_COUNT] = {6, 12, 18, 500, 1000};

double bisectionResultSegments[SEGMENT_COUNT][2] = {{0, 0}, {0, 0}, {0, 0}};

int bisectionTotalSteps = 0;

double discreteNewtonResults[SEGMENT_COUNT] = {0, 0, 0};

int discreteNewtonFailures = 0;

int discreteNewtonTotalSteps[SEGMENT_COUNT] = {0, 0, 0};

double newtonResults[SEGMENT_COUNT] = {0, 0, 0};

int newtonFailures = 0;

int newtonTotalSteps[SEGMENT_COUNT] = {0, 0, 0};

double interpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT];

double chebyshevInterpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT];

double cubicInterpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT];

double Function (double x)
{
    return (pow (x, 9) + M_PI) * cos (log (pow (x, 2) + 1)) / exp (pow (x, 2)) - x / 2019;
}

double FunctionDerivative (double x)
{
    return -(2 * pow (M_E, -pow (x, 2)) * (pow (x, 9) + M_PI) * x *
        sin (log (pow (x, 2) + 1))) / (pow (x, 2) + 1) - 2 *
        pow (M_E, -pow (x, 2)) * (pow (x, 9) + M_PI) * x * cos (log (pow (x, 2) + 1)) + 9 *
        pow (M_E, -pow (x, 2)) * pow (x, 8) * cos (log (pow (x, 2) + 1)) - 1.0 / 2019;
}

void DoBisections (double segments[SEGMENT_COUNT][2])
{
    printf ("Bisections solver.\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        printf ("    Segment index: %d.\n", index);
        int iterations;
        Bisection (Function, &(segments[index][0]), &(segments[index][1]), &iterations);

        bisectionTotalSteps += iterations;
        bisectionResultSegments[index][0] = segments[index][0];
        bisectionResultSegments[index][1] = segments[index][1];

        printf ("        Bisection result segment: [%20.16lf; %20.16lf],\n", segments[index][0], segments[index][1]);
        printf ("        Bisection steps: %d.\n", iterations);
    }
}

void DoDiscreteNewton (double segments[SEGMENT_COUNT][2])
{
    printf ("Discrete Newton solver.\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        printf ("    Segment index: %d.\n", index);
        double x;
        int iterations;

        if (DiscreteNewton (Function, &x, segments[index][0], segments[index][1], &iterations))
        {
            discreteNewtonTotalSteps[index] += iterations;
            discreteNewtonResults[index] = x;

            printf ("        Discrete newton result: %20.16lf.\n", x);
            printf ("        Discrete newton steps: %d.\n", iterations);
        }
        else
        {
            printf ("        Discrete newton failed!");
            discreteNewtonFailures++;
        }
    }
}

void DoNewton (double segments[SEGMENT_COUNT][2])
{
    printf ("Newton solver.\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        printf ("    Segment index: %d.\n", index);
        double x = segments[index][1];
        int iterations;

        if (Newton (Function, FunctionDerivative, &x, segments[index][0], segments[index][1], &iterations))
        {
            newtonTotalSteps[index] += iterations;
            newtonResults[index] = x;

            printf ("        Newton result: %20.16lf.\n", x);
            printf ("        Newton steps: %d.\n", iterations);
        }
        else
        {
            printf ("        Newton failed!\n");
            newtonFailures++;
        }
    }
}

static void PrintInterpolationCoefficients (int nodeSetIndex, int runIndex, double *result)
{
    if (nodeSetIndex <= INTERPOLATION_RESULT_MUTE_AFTER)
    {
        printf ("    Result coefficients for %d nodes, run #%d.\n", interpolationNodeSet[nodeSetIndex], runIndex);
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex]; ++index)
        {
            printf ("        # a%d = %lf\n", index, result[index]);
        }
    }
    else if (runIndex == 0)
    {
        printf ("    Output muted for %d nodes test.\n", interpolationNodeSet[nodeSetIndex]);
    }
}

void DoDifferentialInterpolation ()
{
    printf ("Differential interpolation.\n");
    for (int nodeSetIndex = 0; nodeSetIndex < INTERPOLATION_NODE_SET_COUNT; ++nodeSetIndex)
    {
        time_t totalTime = 0;
        const int runCount = 100;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double *result;
            time_t begin = clock ();
            DifferentialInterpolation (Function, -4, 4, NULL, interpolationNodeSet[nodeSetIndex], &result);

            totalTime += clock () - begin;
            PrintInterpolationCoefficients (nodeSetIndex, runIndex, result);
            free (result);
        }

        interpolationAverageTimeMs[nodeSetIndex] = ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
    }
}

void DoChebyshevDifferentialInterpolation ()
{
    printf ("Chebyshev Differential interpolation.\n");
    for (int nodeSetIndex = 0; nodeSetIndex < INTERPOLATION_NODE_SET_COUNT; ++nodeSetIndex)
    {
        time_t totalTime = 0;
        const int runCount = 100;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double *result;
            time_t begin = clock ();
            ChebyshevDifferentialInterpolation (Function, -4, 4, interpolationNodeSet[nodeSetIndex], &result);

            totalTime += clock () - begin;
            PrintInterpolationCoefficients (nodeSetIndex, runIndex, result);
            free (result);
        }

        chebyshevInterpolationAverageTimeMs[nodeSetIndex] = ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
    }
}

static void PrintSplineCoefficients (int nodeSetIndex, int runIndex, double *result)
{
    if (nodeSetIndex <= INTERPOLATION_RESULT_MUTE_AFTER)
    {
        printf ("    Python splines for %d nodes, run #%d.\n", interpolationNodeSet[nodeSetIndex], runIndex);
        printf ("left = -4\n"
                "right = 4\n"
                "part = (right - left) / %d\n", interpolationNodeSet[nodeSetIndex] - 1);
        printf("nodes = [left + part * i for i in range(%d)]\n", interpolationNodeSet[nodeSetIndex]);

        printf ("a = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%lf", result[4 * index]);
        }

        printf ("]\nb = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%lf", result[4 * index + 1]);
        }

        printf ("]\nc = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%lf", result[4 * index + 2]);
        }

        printf ("]\nd = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%lf", result[4 * index + 3]);
        }

        printf ("]\nxs = generate_xs(left, right, 1000)\n"
                "spline = build_spline(a, b, c, d)\n"
                "ys = [spline(nodes, x) for x in xs]\n"
                "show_plot(xs,ys,show_f=True)\n");
    }
    else if (runIndex == 0)
    {
        printf ("    Output muted for %d nodes test.\n", interpolationNodeSet[nodeSetIndex]);
    }
}

void DoCubicSpline ()
{
    printf ("Cubic Spline\n");
    for (int nodeSetIndex = 0; nodeSetIndex < INTERPOLATION_NODE_SET_COUNT; ++nodeSetIndex)
    {
        time_t totalTime = 0;
        const int runCount = 1;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double *result;
            time_t begin = clock ();
            CubicSpline (Function, -4, 4, interpolationNodeSet[nodeSetIndex] - 1, &result);

            totalTime += clock () - begin;
            PrintSplineCoefficients (nodeSetIndex, runIndex, result);
            free (result);
        }

        cubicInterpolationAverageTimeMs[nodeSetIndex] = ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
    }
}

void PrintReport (FILE *output)
{
    fprintf (output, "#1\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        fprintf (output, "    Bisection result segment %d: [%20.16f; %20.16lf].\n", index,
                 bisectionResultSegments[index][0], bisectionResultSegments[index][1]);
    }

    fprintf (output, "    Bisection steps: %d.\n", bisectionTotalSteps / SEGMENT_COUNT);

    fprintf (output, "#2\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        fprintf (output, "    Discrete newton result segment %d: %22.16lf.\n", index,
                 discreteNewtonResults[index]);
        fprintf (output, "    Discrete newton steps: %d.\n", discreteNewtonTotalSteps[index]);
    }

    fprintf (output, "    Discrete newton failures: %d.\n", discreteNewtonFailures);

    fprintf (output, "#3\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        fprintf (output, "    Newton result segment %d: %22.16lf.\n", index, newtonResults[index]);
        fprintf (output, "    Newton steps: %d.\n", newtonTotalSteps[index]);
    }

    fprintf (output, "    Newton failures: %d.\n", newtonFailures);

    fprintf (output, "#4\n");
    for (int index = 0; index < INTERPOLATION_NODE_SET_COUNT; ++index)
    {
        fprintf (output, "    Average interpolation time for %d nodes: %lfms.\n",
                 interpolationNodeSet[index], interpolationAverageTimeMs[index]);
    }

    fprintf (output, "#5\n");
    for (int index = 0; index < INTERPOLATION_NODE_SET_COUNT; ++index)
    {
        fprintf (output, "    Average Chebyshev interpolation time for %d nodes: %lfms.\n",
                 interpolationNodeSet[index], chebyshevInterpolationAverageTimeMs[index]);
    }

    fprintf (output, "#6\n");
    for (int index = 0; index < INTERPOLATION_NODE_SET_COUNT; ++index)
    {
        fprintf (output, "    Average cubic interpolation time for %d nodes: %lfms.\n",
                 interpolationNodeSet[index], cubicInterpolationAverageTimeMs[index]);
    }
}

int main ()
{
    double segments[SEGMENT_COUNT][2] = {{-2, -1.6}, {-1.4, -1}, {1.7, 2.1}};
    DoBisections (segments);
    DoDiscreteNewton (segments);
    DoNewton (segments);

    printf ("Python support code:\n%s\n", pythonSupportCode);

    DoDifferentialInterpolation ();
    DoChebyshevDifferentialInterpolation ();
    DoCubicSpline ();

    PrintReport (stdout);
    FILE *report = fopen ("report.txt", "w");
    PrintReport (report);
    fclose (report);
    return 0;
}
