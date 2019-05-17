#ifndef __MTWISTER_H
#define __MTWISTER_H

#define STATE_VECTOR_LENGTH 624
#define STATE_VECTOR_M      397 /* changes to STATE_VECTOR_LENGTH also require changes to this */

typedef struct tagMTRand {
    unsigned long mt[STATE_VECTOR_LENGTH];
    int index;
} MTRand;

MTRand SeedRand (unsigned long seed);
unsigned long GenRandLong (MTRand *rand);
double GenRand (MTRand *rand);

#endif /* #ifndef __MTWISTER_H */
