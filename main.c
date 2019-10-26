#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "bisection.h"
#include "discretenewton.h"
#include "newton.h"
#include "differentialinterpolation.h"

#define SEGMENT_COUNT 3

#define INTERPOLATION_NODE_SET_COUNT 6

#define INTERPOLATION_RESULT_MUTE_AFTER 2

int interpolationNodeSet[INTERPOLATION_NODE_SET_COUNT] = {6, 12, 18, 250, 500, 1000};

double bisectionResultSegments[SEGMENT_COUNT][2] = {{0, 0}, {0, 0}, {0, 0}};

int bisectionTotalSteps = 0;

double discreteNewtonResults[SEGMENT_COUNT] = {0, 0, 0};

int discreteNewtonFailures = 0;

int discreteNewtonTotalSteps[SEGMENT_COUNT] = {0, 0, 0};

double newtonResults[SEGMENT_COUNT] = {0, 0, 0};

int newtonFailures = 0;

int newtonTotalSteps[SEGMENT_COUNT] = {0, 0, 0};

double interpolationAverageTimeMs[INTERPOLATION_NODE_SET_COUNT];

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

void DoDifferentialInterpolation()
{
    printf ("Differential interpolation.\n");
    for (int nodeSetIndex = 0; nodeSetIndex < INTERPOLATION_NODE_SET_COUNT; ++nodeSetIndex)
    {
        time_t totalTime = 0;
        const int runCount = 1000;

        for (int runIndex = 0; runIndex < runCount; ++runIndex)
        {
            double *result;
            time_t begin = clock ();
            DifferentialInterpolation (Function, -4, 4, interpolationNodeSet[nodeSetIndex], &result);
            totalTime += clock () - begin;

            if (nodeSetIndex <= INTERPOLATION_RESULT_MUTE_AFTER)
            {
                printf ("    Result coefficients for %d nodes, run #%d.\n",
                        interpolationNodeSet[nodeSetIndex], runIndex);

                for (int index = 0; index < interpolationNodeSet[nodeSetIndex]; ++index)
                {
                    printf ("        # a%d = %lf\n", index, result[index]);
                }
            }
            else if (runIndex == 0)
            {
                printf ("Output muted for %d nodes test.\n", interpolationNodeSet[nodeSetIndex]);
            }

            free (result);
        }

        interpolationAverageTimeMs[nodeSetIndex] = ((double) totalTime / runCount) / (CLOCKS_PER_SEC / 1000.0);
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
        fprintf (output, "    Average time for %d nodes: %lfms.\n",
            interpolationNodeSet[index], interpolationAverageTimeMs[index]);
    }
}

int main ()
{
    double segments[SEGMENT_COUNT][2] = {{-2, -1.6}, {-1.4, -1}, {1.7, 2.1}};
    DoBisections (segments);
    DoDiscreteNewton (segments);
    DoNewton (segments);
    DoDifferentialInterpolation ();

    PrintReport (stdout);
    FILE *report = fopen ("report.txt", "w");
    PrintReport (report);
    fclose (report);
    return 0;
}
