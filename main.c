#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

#include "matrixutils.h"
#include "bisection.h"
#include "discretenewton.h"
#include "newton.h"
#include "differentialinterpolation.h"
#include "cubicspline.h"
#include "bezier.h"
#include "quadraticregression.h"
#include "bilagrange.h"
#include "simpleintegrals.h"
#include "gaussianintegral.h"
#include "runge.h"
#include "bicubicspline.h"

#define SEGMENT_COUNT 3

#define INTERPOLATION_NODE_SET_COUNT 5

#define INTERPOLATION_RESULT_MUTE_AFTER 2

#define REGRESSION_POWER_SET_COUNT 9

#define I 4.967532679086564

// FIXME: Doing now only up to 9, must do up to 21!
#define INTEGRAL_STEP_POWERS 10

#define NEWTON_POWERS_COUNT 6

#define GAUSSIAN_POWERS_COUNT 6

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
    "    return sorted([rnd.uniform(left,right) for i in range(n)])\n"
    "\n"
    "xleft, xright, yleft, yright = -4, 4, -6, 6\n"
    "def build_mesh(zs):\n"
    "    x_step = (xright - xleft)/math.ceil(math.sqrt(len(zs)))\n"
    "    y_step = (yright - yleft)/math.ceil(math.sqrt(len(zs)))  \n"
    "    plt.style.use('ggplot')\n"
    "    fig = plt.figure()\n"
    "    ax = Axes3D(fig)\n"
    "    ax.view_init(30, 30)\n"
    "    ax.set_zlim(-3, 3)\n"
    "    plt.rcParams['figure.figsize'] = (15, 5)\n"
    "    Xs = [xleft + x_step * i for i in range(math.ceil(math.sqrt(len(zs))))]\n"
    "    Ys = [yleft + y_step * i for i in range(math.ceil(math.sqrt(len(zs))))]\n"
    "    X, Y = np.meshgrid(Xs, Ys)\n"
    "    \n"
    "    Z = []\n"
    "    for row in range(0, math.ceil(math.sqrt(len(zs)))):\n"
    "        R = []\n"
    "        for col in range(0, math.ceil(math.sqrt(len(zs)))):\n"
    "            R.append(zs[col * math.ceil(math.sqrt(len(zs))) + row])\n"
    "        Z.append(R)\n"
    "    \n"
    "    ax.plot_surface(X, Y, np.array(Z))\n"
    "    plt.grid(True)\n"
    "    plt.show()\n"
    "\n";

const int interpolationNodeSet[INTERPOLATION_NODE_SET_COUNT] = {6, 12, 18, 500, 1000};

const int regressionPowerSet[REGRESSION_POWER_SET_COUNT] = {1, 2, 3, 4, 5, 6, 9, 12, 15};

const int newtonPowers[NEWTON_POWERS_COUNT] = {5, 7, 9, 11, 13, 15};

const int gaussianPowers[GAUSSIAN_POWERS_COUNT] = {2, 3, 4, 5, 6, 7};

double bisectionResultSegments[SEGMENT_COUNT][2] = {{0, 0}, {0, 0}, {0, 0}};

int bisectionTotalSteps = 0;

double discreteNewtonResults[SEGMENT_COUNT] = {0, 0, 0};

int discreteNewtonFailures = 0;

int discreteNewtonTotalSteps[SEGMENT_COUNT] = {0, 0, 0};

double newtonResults[SEGMENT_COUNT] = {0, 0, 0};

int newtonFailures = 0;

int newtonTotalSteps[SEGMENT_COUNT] = {0, 0, 0};

double interpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT] = {0};

double chebyshevInterpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT] = {0};

double cubicInterpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT] = {0};

double bezierInterpolationAverageTimeMs = 0;

double quadraticRegressionInterpolationAverageTimeMs[REGRESSION_POWER_SET_COUNT] = {0};

double biLagrangeInterpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT] = {0};

double bicubicInterpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT] = {0};

double leftRectangleI[INTEGRAL_STEP_POWERS] = {0};

double rightRectangleI[INTEGRAL_STEP_POWERS] = {0};

double mediumRectangleI[INTEGRAL_STEP_POWERS] = {0};

double trapeziumI[INTEGRAL_STEP_POWERS] = {0};

double simpsonI[INTEGRAL_STEP_POWERS] = {0};

double newtonI[INTEGRAL_STEP_POWERS * NEWTON_POWERS_COUNT] = {0};

double gaussianI[INTEGRAL_STEP_POWERS * GAUSSIAN_POWERS_COUNT] = {0};

double rungeNewtonI[NEWTON_POWERS_COUNT] = {0};

double rungeGaussianI[GAUSSIAN_POWERS_COUNT] = {0};

double leftRectangleAverageMs[INTEGRAL_STEP_POWERS] = {0};

double rightRectangleAverageMs[INTEGRAL_STEP_POWERS] = {0};

double mediumRectangleAverageMs[INTEGRAL_STEP_POWERS] = {0};

double trapeziumAverageMs[INTEGRAL_STEP_POWERS] = {0};

double simpsonAverageMs[INTEGRAL_STEP_POWERS] = {0};

double newtonAverageMs[INTEGRAL_STEP_POWERS * NEWTON_POWERS_COUNT] = {0};

double gaussianAverageMs[INTEGRAL_STEP_POWERS * GAUSSIAN_POWERS_COUNT] = {0};

double rungeNewtonAverageMs[INTEGRAL_STEP_POWERS * NEWTON_POWERS_COUNT] = {0};

double rungeGaussianAverageMs[INTEGRAL_STEP_POWERS * GAUSSIAN_POWERS_COUNT] = {0};

unsigned long rungeNewtonPartitions[NEWTON_POWERS_COUNT] = {0};

unsigned long rungeNewtonSteps[NEWTON_POWERS_COUNT] = {0};

unsigned long rungeGaussianPartitions[GAUSSIAN_POWERS_COUNT] = {0};

unsigned long rungeGaussianSteps[GAUSSIAN_POWERS_COUNT] = {0};

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

double GFunction (double x, double y)
{
    return (pow (x, 9) + M_PI) * cos (log (pow (y, 4) + 1)) / ((pow (y, 2) + M_E) * exp (x * x));
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

static void PrintInterpolationCoefficients (int nodeSetIndex, int runIndex, double *result, double *xs)
{
    if (nodeSetIndex <= INTERPOLATION_RESULT_MUTE_AFTER)
    {
        if (runIndex != 0)
        {
            return;
        }

        printf ("    Result formula for %d nodes.\n", interpolationNodeSet[nodeSetIndex]);
        double step = 8.0 / (interpolationNodeSet[nodeSetIndex] - 1);

        for (int index = 0; index < interpolationNodeSet[nodeSetIndex]; ++index)
        {
            printf ("%0.16lf", result[index]);
            if (index > 0)
            {
                for (int mulIndex = 0; mulIndex < index; ++mulIndex)
                {
                    if (xs == NULL)
                    {
                        printf ("*(x-%20.16lf)", -4 + step * mulIndex);
                    }
                    else
                    {
                        printf ("*(x-%20.16lf)", xs[mulIndex]);
                    }
                }
            }

            if (index < interpolationNodeSet[nodeSetIndex] - 1)
            {
                printf ("+");
            }
            else
            {
                printf ("\n");
            }
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
        const int runCount = 20;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double *result;
            time_t begin = clock ();
            DifferentialInterpolation (Function, -4, 4, NULL, interpolationNodeSet[nodeSetIndex], &result);

            totalTime += clock () - begin;
            PrintInterpolationCoefficients (nodeSetIndex, runIndex, result, NULL);
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
        const int runCount = 20;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double *result;
            time_t begin = clock ();
            ChebyshevDifferentialInterpolation (Function, -4, 4, interpolationNodeSet[nodeSetIndex], &result);

            totalTime += clock () - begin;
            double *xs = calloc (interpolationNodeSet[nodeSetIndex], sizeof (double));

            for (int index = 0; index < interpolationNodeSet[nodeSetIndex]; ++index)
            {
                xs[index] = cos ((2 * index + 1) * M_PI / (2 * interpolationNodeSet[nodeSetIndex] + 2)) * 4;
            }

            PrintInterpolationCoefficients (nodeSetIndex, runIndex, result, xs);
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
        printf ("nodes = [left + part * i for i in range(%d)]\n", interpolationNodeSet[nodeSetIndex]);

        printf ("a = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%0.16lf", result[4 * index]);
        }

        printf ("]\nb = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%0.16lf", result[4 * index + 1]);
        }

        printf ("]\nc = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%0.16lf", result[4 * index + 2]);
        }

        printf ("]\nd = [");
        for (int index = 0; index < interpolationNodeSet[nodeSetIndex] - 1; ++index)
        {
            if (index != 0)
            {
                printf (", ");
            }

            printf ("%0.16lf", result[4 * index + 3]);
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

static void DoCubicSpline ()
{
    printf ("Cubic Spline.\n");
    for (int nodeSetIndex = 0; nodeSetIndex < INTERPOLATION_NODE_SET_COUNT; ++nodeSetIndex)
    {
        time_t totalTime = 0;
        const int runCount = 10;

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

static void PrintPythonArray (const char *name, double *values, int size)
{
    printf ("%s = [", name);
    for (int index = 0; index < size; ++index)
    {
        if (index != 0)
        {
            printf (", ");
        }

        printf ("%0.16lf", values[index]);
    }

    printf ("]\n");
}

static void PrintBezierAsTFunction (double *bezier, int count)
{
    for (int index = 0; index < count; ++index)
    {
        if (index > 0)
        {
            printf (" + ");
        }

        printf ("%0.16lf", bezier[index]);
        if (index > 1)
        {
            printf ("*(t**%d)", index);
        }
        else if (index == 1)
        {
            printf ("*t");
        }
    }
}

void PrintBezier (double *x, double *y, int count, double *resultX, double *resultY)
{
    printf ("    Bezier curve for %d nodes python code:\n", count);
    PrintPythonArray ("control_xs", x, count);
    PrintPythonArray ("control_points", y, count);

    printf ("x_func = lambda t: ");
    PrintBezierAsTFunction (resultX, count);
    printf ("\n");

    printf ("y_func = lambda t: ");
    PrintBezierAsTFunction (resultY, count);
    printf ("\n");

    printf ("t_array = generate_xs(0,1,300)\n"
            "bezier_xs = [x_func(t) for t in t_array]\n"
            "bezier_ys = [y_func(t) for t in t_array]\n"
            "plt.scatter(control_xs, control_points)\n");
    printf ("for i in range(%d):\n"
            "    plt.annotate(i, (control_xs[i], control_points[i]))\n", count);
    printf ("plt.plot(bezier_xs, bezier_ys)\n");
}

void DoBezier ()
{
    printf ("Bezier.\n");
    time_t totalTime = 0;
    const int runCount = 30;

    for (int runIndex = 0; runIndex < runCount; ++runIndex)
    {
        const int count = 40;
        double *x = calloc (count, sizeof (double));
        double *y = calloc (count, sizeof (double));

        for (int index = 0; index < count; ++index)
        {
            x[index] = -4 + rand () * 8.0 / RAND_MAX;
            // Если нужна нормальная кривая, то берём отсортированные точки.
            x[index] = -4 + 8.0 / (count - 1) * index;
            y[index] = Function (x[index]);
        }

        double *resultX;
        double *resultY;

        time_t begin = clock ();
        BezierCasteljau (x, y, count, &resultX, &resultY);
        totalTime += clock () - begin;
        PrintBezier (x, y, count, resultX, resultY);

        free (x);
        free (y);
        free (resultX);
        free (resultY);
    }

    bezierInterpolationAverageTimeMs = ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
}

static void PrintQuadraticRegressionAsXFunction (double **quadratic, int count)
{
    for (int index = 0; index < count; ++index)
    {
        if (index > 0)
        {
            printf (" + ");
        }

        printf ("%0.16lf", quadratic[index][0]);
        if (index > 1)
        {
            printf ("*(x**%d)", index);
        }
        else if (index == 1)
        {
            printf ("*x");
        }
    }
}

void PrintQuadraticRegression (double *x, double *y, int count, double **result, int power, int runIndex)
{
    if (runIndex != 0)
    {
        return;
    }

    printf ("    Quadratic Regression of power %d for %d nodes python code:\n", power, count);
    PrintPythonArray ("control_xs", x, count);
    PrintPythonArray ("control_points", y, count);

    printf ("func = lambda x: ");
    PrintQuadraticRegressionAsXFunction (result, power + 1);
    printf ("\n");

    printf ("x_array = generate_xs(left, right,1000)\n"
            "ys = [func(x) for x in x_array]\n"
            "plt.scatter(control_xs, control_points)\n");
    printf ("plt.plot(x_array, ys)\n");
}

static void DoQuadraticRegression ()
{
    printf ("Quadratic regression.\n");
    const int runCount = 10;

    for (int powerIndex = 0; powerIndex < REGRESSION_POWER_SET_COUNT; ++powerIndex)
    {
        const int power = regressionPowerSet[powerIndex];
        time_t totalTime = 0;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            const int count = 100;
            double *x = calloc (count, sizeof (double));
            double *y = calloc (count, sizeof (double));

            for (int index = 0; index < count; ++index)
            {
                // TODO: Use equidistant points by now (despite lab requirements).
                x[index] = -4 + 8.0 / (count - 1) * index;
                y[index] = Function (x[index]);
            }

            double **result;
            time_t begin = clock ();
            QuadraticRegression (x, y, count, power, &result);
            totalTime += clock () - begin;
            PrintQuadraticRegression (x, y, count, result, power, runIndex);

            free (x);
            free (y);
            FreeMatrix (result, power + 1, 1);
        }

        quadraticRegressionInterpolationAverageTimeMs[powerIndex] =
            ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
    }
}

static void DoBiLagrange ()
{
    printf ("BiLagrange.\n");
    for (int nodeSetIndex = 0; nodeSetIndex <= INTERPOLATION_RESULT_MUTE_AFTER; ++nodeSetIndex)
    {
        time_t totalTime = 0;
        const int runCount = 5;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double *x = calloc (interpolationNodeSet[nodeSetIndex], sizeof (double));
            x[0] = -4;

            for (int index = 1; index < interpolationNodeSet[nodeSetIndex]; ++index)
            {
                x[index] = x[index - 1] + 8.0 / (interpolationNodeSet[nodeSetIndex] - 1);
            }

            double *y = calloc (interpolationNodeSet[nodeSetIndex], sizeof (double));
            y[0] = -6;

            for (int index = 1; index < interpolationNodeSet[nodeSetIndex]; ++index)
            {
                y[index] = y[index - 1] + 12.0 / (interpolationNodeSet[nodeSetIndex] - 1);
            }

            double **z = AllocateMatrix (interpolationNodeSet[nodeSetIndex], interpolationNodeSet[nodeSetIndex]);
            for (int row = 0; row < interpolationNodeSet[nodeSetIndex]; ++row)
            {
                for (int col = 0; col < interpolationNodeSet[nodeSetIndex]; ++col)
                {
                    z[row][col] = GFunction (x[row], y[col]);
                }
            }

            printf ("    Polynomial for %dx%d split: ",
                    interpolationNodeSet[nodeSetIndex], interpolationNodeSet[nodeSetIndex]);
            time_t begin = clock ();

            PrintBiLagrangePolynomial (x, y, interpolationNodeSet[nodeSetIndex], z, stdout);
            totalTime += clock () - begin;
            printf ("\n");

            free (x);
            free (y);
            FreeMatrix (z, interpolationNodeSet[nodeSetIndex], interpolationNodeSet[nodeSetIndex]);
        }

        biLagrangeInterpolationAverageTimeMs[nodeSetIndex] =
            ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
    }
}

void DoIntegrals ()
{
    const int runCount = 10;
    long parts = 1;

    for (int i = 0; i < INTEGRAL_STEP_POWERS; ++i)
    {
        printf ("Integral parts: %ld\n", parts);
        for (int run = 0; run < runCount; ++run)
        {
            time_t begin = clock ();
            leftRectangleI[i] = LeftRectangleIntegral (Function, -4, 4, parts);
            leftRectangleAverageMs[i] += (clock () - begin) * 1.0 / runCount;
        }

        for (int run = 0; run < runCount; ++run)
        {
            time_t begin = clock ();
            rightRectangleI[i] = RightRectangleIntegral (Function, -4, 4, parts);
            rightRectangleAverageMs[i] += (clock () - begin) * 1.0 / runCount;
        }

        for (int run = 0; run < runCount; ++run)
        {
            time_t begin = clock ();
            mediumRectangleI[i] = MediumRectangleIntegral (Function, -4, 4, parts);
            mediumRectangleAverageMs[i] += (clock () - begin) * 1.0 / runCount;
        }

        for (int run = 0; run < runCount; ++run)
        {
            time_t begin = clock ();
            trapeziumI[i] = TrapeziumRectangleIntegral (Function, -4, 4, parts);
            trapeziumAverageMs[i] += (clock () - begin) * 1.0 / runCount;
        }

        for (int run = 0; run < runCount; ++run)
        {
            time_t begin = clock ();
            simpsonI[i] = SimpsonRectangleIntegral (Function, -4, 4, parts);
            simpsonAverageMs[i] += (clock () - begin) * 1.0 / runCount;
        }

        for (int newtonPowerIndex = 0; newtonPowerIndex < NEWTON_POWERS_COUNT; ++newtonPowerIndex)
        {
            int index = i * NEWTON_POWERS_COUNT + newtonPowerIndex;
            time_t begin = clock ();
            newtonI[index] = NewtonIntegral (Function, newtonPowers[newtonPowerIndex], -4, 4, parts);
            newtonAverageMs[index] += (clock () - begin) * 1.0 / runCount;
        }

        for (int gaussianPowerIndex = 0; gaussianPowerIndex < GAUSSIAN_POWERS_COUNT; ++gaussianPowerIndex)
        {
            int index = i * GAUSSIAN_POWERS_COUNT + gaussianPowerIndex;
            time_t begin = clock ();
            gaussianI[index] = GaussianIntegral (Function, gaussianPowers[gaussianPowerIndex], -4, 4, parts);
            gaussianAverageMs[index] += (clock () - begin) * 1.0 / runCount;
        }

        parts *= 4;
    }
}

void DoRungeRule ()
{
    printf ("Running with Runge rule...\n");
    const int runCount = 10;
    const double epsilon = pow (10, -14);

    for (int run = 0; run < runCount; ++run)
    {
        printf ("    Iteration %d...\n", run);
        for (int newtonPowerIndex = 0; newtonPowerIndex < NEWTON_POWERS_COUNT; ++newtonPowerIndex)
        {
            time_t begin = clock ();
            rungeNewtonI[newtonPowerIndex] = RungeRule (NewtonIntegral, Function,
                                                        newtonPowers[newtonPowerIndex],
                                                        newtonPowers[newtonPowerIndex], -4, 4, epsilon,
                                                        rungeNewtonPartitions + newtonPowerIndex,
                                                        rungeNewtonSteps + newtonPowerIndex);
            rungeNewtonAverageMs[newtonPowerIndex] += (clock () - begin) * 1.0 / runCount;
        }

        for (int gaussianPowerIndex = 0; gaussianPowerIndex < GAUSSIAN_POWERS_COUNT; ++gaussianPowerIndex)
        {
            time_t begin = clock ();
            rungeGaussianI[gaussianPowerIndex] = RungeRule (GaussianIntegral, Function,
                                                            gaussianPowers[gaussianPowerIndex],
                                                            gaussianPowers[gaussianPowerIndex] * 2, -4, 4,  epsilon,
                                                            rungeGaussianPartitions + gaussianPowerIndex,
                                                            rungeGaussianSteps + gaussianPowerIndex);
            rungeGaussianAverageMs[gaussianPowerIndex] += (clock () - begin) * 1.0 / runCount;
        }
    }
}

void DoBiCubic ()
{
    printf ("BiCubic.\n");
    const int gridSize = 250;
    double *Zs = calloc (gridSize * gridSize, sizeof (double));
    double stepX = 8.0 / (gridSize - 1);
    double stepY = 12.0 / (gridSize - 1);

    double x = -4;
    double y = -6;

    for (int row = 0; row < gridSize; ++row)
    {
        for (int col = 0; col < gridSize; ++col)
        {
            Zs[row * gridSize + col] = GFunction (x, y);
            y += stepY;
        }

        y = -6;
        x += stepX;
    }

    printf ("    Points of real function: ");
    PrintPythonArray ("Zs", Zs, gridSize * gridSize);
    printf ("\n");

    for (int nodeSetIndex = 0; nodeSetIndex < INTERPOLATION_NODE_SET_COUNT; ++nodeSetIndex)
    {
        time_t totalTime = 0;
        const int runCount = 10;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double **xRes;
            time_t begin = clock ();
            BiCubicSpline (GFunction, -4, 4, -6, 6, interpolationNodeSet[nodeSetIndex] + 1, &xRes);
            totalTime += clock () - begin;

            x = -4;
            y = -6;

            for (int row = 0; row < gridSize; ++row)
            {
                for (int col = 0; col < gridSize; ++col)
                {
                    Zs[row * gridSize + col] = BiCubicSplineCalculate(
                        x, y, -4, 4, -6, 6, interpolationNodeSet[nodeSetIndex] + 1, xRes);
                    y += stepY;
                }

                y = -6;
                x += stepX;
            }

            if (runIndex == 0)
            {
                printf ("    Points for %dx%d split mesh: ",
                        interpolationNodeSet[nodeSetIndex], interpolationNodeSet[nodeSetIndex]);
                PrintPythonArray ("Zs", Zs, gridSize * gridSize);
                printf ("\n");
            }

            for (int index = 0; index < interpolationNodeSet[nodeSetIndex] * 4; ++index)
            {
                free(xRes[index]);
            }

            free (xRes);
        }

        bicubicInterpolationAverageTimeMs[nodeSetIndex] =
            ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
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

    fprintf (output, "#7\n");
    fprintf (output, "    Bezier interpolation time for 40 nodes: %lfms.\n", bezierInterpolationAverageTimeMs);

    fprintf (output, "#8\n");
    for (int index = 0; index < REGRESSION_POWER_SET_COUNT; ++index)
    {
        fprintf (output, "    Average quadratic regression time for power %d: %lfms.\n",
                 regressionPowerSet[index], quadraticRegressionInterpolationAverageTimeMs[index]);
    }

    fprintf (output, "#9\n");
    for (int index = 0; index <= INTERPOLATION_RESULT_MUTE_AFTER; ++index)
    {
        fprintf (output, "    Average bilagrange interpolation time for %d nodes: %lfms.\n",
                 interpolationNodeSet[index], biLagrangeInterpolationAverageTimeMs[index]);
    }

    fprintf (output, "#10\n");
    for (int index = 0; index < INTERPOLATION_NODE_SET_COUNT; ++index)
    {
        fprintf (output, "    Average bicubicinterpolation time for %d split: %lfms.\n",
                 interpolationNodeSet[index], bicubicInterpolationAverageTimeMs[index]);
    }

    fprintf (output, "#11\n");
    for (int i = 0; i < INTEGRAL_STEP_POWERS; ++i)
    {
        fprintf (output, "    Power: %d.\n", i);
        fprintf (output, "        Left rectangle difference: %22.16lf.\n", fabs (I - leftRectangleI[i]));
        fprintf (output, "        Right rectangle difference: %22.16lf.\n", fabs (I - rightRectangleI[i]));
        fprintf (output, "        Medium rectangle difference: %22.16lf.\n", fabs (I - mediumRectangleI[i]));
        fprintf (output, "        Trapezium difference: %22.16lf.\n", fabs (I - trapeziumI[i]));
        fprintf (output, "        Simpson difference: %22.16lf.\n", fabs (I - simpsonI[i]));

        fprintf (output, "        Left rectangle time: %22.16lfms.\n", leftRectangleAverageMs[i]);
        fprintf (output, "        Right rectangle time: %22.16lfms.\n", rightRectangleAverageMs[i]);
        fprintf (output, "        Medium rectangle time: %22.16lfms.\n", mediumRectangleAverageMs[i]);
        fprintf (output, "        Trapezium time: %22.16lfms.\n", trapeziumAverageMs[i]);
        fprintf (output, "        Simpson time: %22.16lfms.\n", simpsonAverageMs[i]);

        for (int newtonIndex = 0; newtonIndex < NEWTON_POWERS_COUNT; ++newtonIndex)
        {
            int index = i * NEWTON_POWERS_COUNT + newtonIndex;
            fprintf (output, "        Newton power %d difference: %22.16lf.\n",
                     newtonPowers[newtonIndex], fabs (I - newtonI[index]));
            fprintf (output, "        Newton power %d time: %22.16lfms.\n",
                     newtonPowers[newtonIndex], newtonAverageMs[index]);
        }

        for (int gaussianIndex = 0; gaussianIndex < GAUSSIAN_POWERS_COUNT; ++gaussianIndex)
        {
            int index = i * GAUSSIAN_POWERS_COUNT + gaussianIndex;
            fprintf (output, "        Gaussian power %d difference: %22.16lf.\n",
                     gaussianPowers[gaussianIndex], fabs (I - gaussianI[index]));
            fprintf (output, "        Gaussian power %d time: %22.16lfms.\n",
                     gaussianPowers[gaussianIndex], gaussianAverageMs[index]);
        }
    }

    fprintf (output, "#12\n");
    for (int newtonIndex = 0; newtonIndex < NEWTON_POWERS_COUNT; ++newtonIndex)
    {
        fprintf (output, "    Runge Newton power %d difference: %22.16lf.\n",
                 newtonPowers[newtonIndex], fabs (I - rungeNewtonI[newtonIndex]));
        fprintf (output, "    Runge Newton power %d time: %22.16lfms.\n",
                 newtonPowers[newtonIndex], rungeNewtonAverageMs[newtonIndex]);
        fprintf (output, "    Runge Newton power %d partitions: %ld.\n",
                 newtonPowers[newtonIndex], rungeNewtonPartitions[newtonIndex]);
        fprintf (output, "    Runge Newton power %d steps: %ld.\n",
                 newtonPowers[newtonIndex], rungeNewtonSteps[newtonIndex]);
    }

    for (int gaussianIndex = 0; gaussianIndex < GAUSSIAN_POWERS_COUNT; ++gaussianIndex)
    {
        fprintf (output, "    Runge Gaussian power %d difference: %22.16lf.\n",
                 gaussianPowers[gaussianIndex], fabs (I - rungeGaussianI[gaussianIndex]));
        fprintf (output, "    Runge Gaussian power %d time: %22.16lfms.\n",
                 gaussianPowers[gaussianIndex], rungeGaussianAverageMs[gaussianIndex]);
        fprintf (output, "    Runge Gaussian power %d partitions: %ld.\n",
                 gaussianPowers[gaussianIndex], rungeGaussianPartitions[gaussianIndex]);
        fprintf (output, "    Runge Gaussian power %d steps: %ld.\n",
                 gaussianPowers[gaussianIndex], rungeGaussianSteps[gaussianIndex]);
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
    DoBezier ();
    DoQuadraticRegression ();
    DoBiLagrange ();
    DoIntegrals ();
    DoRungeRule ();
    DoBiCubic ();

    PrintReport (stdout);
    FILE *report = fopen ("report.txt", "w");
    PrintReport (report);
    fclose (report);
    return 0;
}
