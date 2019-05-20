#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "def.h"
#include "mtwister.h"
#include "matrixutils.h"
#include "powervalue.h"
#include "bisection.h"
#include "discretenewton.h"

#define SEGMENT_COUNT 3

extern MTRand *GlobalRand;

double powerMinNorm = INFINITY;

double powerMaxNorm = 0.0;

double powerAverageNorm = 0.0;

clock_t powerTotalTime = 0;

double bisectionResultSegments[SEGMENT_COUNT][2] = {{0, 0}, {0, 0}, {0, 0}};

int bisectionTotalSteps = 0;

double discreteNewtonResultsSum[SEGMENT_COUNT] = {0, 0, 0};

int descreteNewtonFailures = 0;

int discreteNewtonTotalSteps[SEGMENT_COUNT] = {0, 0, 0};

double Function (double x)
{
    return (pow (x, 9) + M_PI) * cos (log (pow (x, 2) + 1)) / exp (pow (x, 2)) - x / 2018;
}

void TestPowerMethod (double **A)
{
    double **myA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **vector;
    double **multiplicationResult = AllocateMatrix (MATRIX_SIZE, 1);

    clock_t begin = clock ();
    double value = PowerValue (myA, MATRIX_SIZE, &vector);
    powerTotalTime += clock () - begin;

    MultiplyMatrices (myA, vector, multiplicationResult, MATRIX_SIZE, MATRIX_SIZE, 1);
    double normal = 0.0;

    for (int index = 0; index < MATRIX_SIZE; ++index)
    {
        normal += pow (vector[index][0] * value - multiplicationResult[index][0], 2);
    }

    normal = sqrt (normal);
    printf ("Power method first value: %20.16lf.\n", value);
    printf ("Power method first value search error: %20.16lf.\n", normal);

    powerMinNorm = fmin (powerMinNorm, normal);
    powerMaxNorm = fmax (powerMaxNorm, normal);
    powerAverageNorm += normal / RUN_COUNT;

    double **myB = CopyMatrix (myA, MATRIX_SIZE, MATRIX_SIZE);
    for (int index = 0; index < MATRIX_SIZE; ++index)
    {
        myB[index][index] -= value;
    }

    begin = clock ();
    double secondValue = PowerValue (myB, MATRIX_SIZE, &vector) + value;
    powerTotalTime += clock () - begin;

    MultiplyMatrices (myA, vector, multiplicationResult, MATRIX_SIZE, MATRIX_SIZE, 1);
    normal = 0.0;

    for (int index = 0; index < MATRIX_SIZE; ++index)
    {
        normal += pow (vector[index][0] * secondValue - multiplicationResult[index][0], 2);
    }

    normal = sqrt (normal);
    printf ("Power method second value: %20.16lf.\n", secondValue);
    printf ("Power method second value search error: %20.16lf.\n", normal);

    powerMinNorm = fmin (powerMinNorm, normal);
    powerMaxNorm = fmax (powerMaxNorm, normal);
    powerAverageNorm += normal / RUN_COUNT;
}

void DoBisections (double segments[SEGMENT_COUNT][2])
{
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        printf ("Segment index: %d.\n", index);
        int iterations;
        Bisection (Function, &(segments[index][0]), &(segments[index][1]), &iterations);

        bisectionTotalSteps += iterations;
        bisectionResultSegments[index][0] += segments[index][0] / RUN_COUNT;
        bisectionResultSegments[index][1] += segments[index][1] / RUN_COUNT;

        printf ("Bisection result segment: [%7.5lf; %7.5lf],\n", segments[index][0], segments[index][1]);
        printf ("Bisection steps: %d.\n", iterations);
    }
}

void DoDiscreteNewton (double segments[SEGMENT_COUNT][2])
{
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        printf ("Segment index: %d.\n", index);
        double x;
        int iterations;

        if (DescriteNewton (Function, &x, segments[index][0], segments[index][1], &iterations))
        {
            discreteNewtonTotalSteps[index] += iterations;
            discreteNewtonResultsSum[index] += x;

            printf ("Discrete newton result: %20.16lf.\n", x);
            printf ("Discrete newton steps: %d.\n", iterations);
        }
        else
        {
            printf ("Discrete newton failed!");
            descreteNewtonFailures++;
        }
    }
}

void MainCycle ()
{
    double **A = AllocateMatrix (MATRIX_SIZE, MATRIX_SIZE);
    FillDefaultMatrix (A, MATRIX_SIZE, MATRIX_SIZE);

    TestPowerMethod (A);

    FreeMatrix (A, MATRIX_SIZE, MATRIX_SIZE);

    double segments[SEGMENT_COUNT][2] = {{-2, -1.6}, {-1.4, -1}, {1.7, 2.1}};
    DoBisections (segments);
    DoDiscreteNewton (segments);
}

void PrintReport (FILE *output)
{
    fprintf (output, "#1\n");
    fprintf (output, "    Power method min error normal: %20.16lf.\n", powerMinNorm);
    fprintf (output, "    Power method max error normal: %20.16lf.\n", powerMaxNorm);
    fprintf (output, "    Power method average error normal: %20.16lf.\n", powerAverageNorm);
    fprintf (output, "    Power method average time: %dms.\n\n",
             (int) (powerTotalTime * CLOCKS_PER_SEC / 1000 / RUN_COUNT / 2));

    fprintf (output, "#2\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        fprintf (output, "    Bisection result segment %d: [%7.5lf; %7.5lf].\n", index,
                 bisectionResultSegments[index][0], bisectionResultSegments[index][1]);
    }

    fprintf (output, "    Bisection average steps: %10.7lf.\n", bisectionTotalSteps * 1.0 / RUN_COUNT / SEGMENT_COUNT);

    fprintf (output, "#3\n");
    for (int index = 0; index < SEGMENT_COUNT; ++index)
    {
        fprintf (output, "    Discrete newton result segment %d: %20.16lf.\n", index,
                 discreteNewtonResultsSum[index] / (RUN_COUNT - descreteNewtonFailures));
        
        fprintf (output, "    Discrete newton average steps: %10.7lf.\n",
                discreteNewtonTotalSteps[index] * 1.0 / RUN_COUNT);
    }

    fprintf (output, "    Discrete newton failures: %d.\n", descreteNewtonFailures);
}

int main ()
{
    GlobalRand = malloc (sizeof (MTRand));
    *GlobalRand = SeedRand (1377);

    for (int index = 0; index < RUN_COUNT; ++index)
    {
        printf ("### Run %d ### \n", index);
        MainCycle ();

        printf ("\n\n");
    }

    PrintReport (stdout);
    FILE *report = fopen ("report.txt", "w");
    PrintReport (report);
    fclose (report);

    free (GlobalRand);
    return 0;
}
