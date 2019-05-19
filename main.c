#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "def.h"
#include "mtwister.h"
#include "matrixutils.h"
#include "powervalue.h"

extern MTRand *GlobalRand;

double powerMinNorm = INFINITY;

double powerMaxNorm = 0.0;

double powerAverageNorm = 0.0;

clock_t powerTotalTime;

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

void MainCycle ()
{
    double **A = AllocateMatrix (MATRIX_SIZE, MATRIX_SIZE);
    FillDefaultMatrix (A, MATRIX_SIZE, MATRIX_SIZE);

    TestPowerMethod (A);

    FreeMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
}

void PrintReport (FILE *output)
{
    fprintf (output, "#1\n");
    fprintf (output, "    Power method min error normal: %20.16lf.\n", powerMinNorm);
    fprintf (output, "    Power method max error normal: %20.16lf.\n", powerMaxNorm);
    fprintf (output, "    Power method average error normal: %20.16lf.\n", powerAverageNorm);
    fprintf (output, "    Power method average time: %dms.\n\n",
             (int) (powerTotalTime * CLOCKS_PER_SEC / 1000 / RUN_COUNT / 2));
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
