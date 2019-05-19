#include "powervalue.h"
#include "matrixutils.h"
#include "def.h"

#include <stdbool.h>
#include <math.h>

#define BARRIER 0.0000000000000001
#define MAX_ITERATION 100000

double PowerValue (double **A, int matrixSize, double ***valueVector)
{
    double **outputVector = AllocateMatrix (matrixSize, 1);
    *valueVector = AllocateMatrix (matrixSize, 1);
    (*valueVector)[0][0] = 1;
    double value = 1.0;
    double previousValues[2] = {0.0, 0.0};
    bool working = true;
    int iteration = 0;

    while (working && iteration < MAX_ITERATION)
    {
        previousValues[0] = previousValues[1];
        previousValues[1] = value;
        MultiplyMatrices (A, *valueVector, outputVector, matrixSize, matrixSize, 1);

        value = 0.0;
        for (int index = 0; index < matrixSize; ++index)
        {
            if (fabs (value) < fabs (outputVector[index][0]))
            {
                value = outputVector[index][0];
            }
        }

        for (int index = 0; index < matrixSize; ++index)
        {
            (*valueVector)[index][0] = outputVector[index][0] / value;
        }

        if (fabs (value - previousValues[0]) < BARRIER)
        {
            working = false;
        }

        iteration++;
    }

    FreeMatrix (outputVector, matrixSize, 1);
    return value;
}
