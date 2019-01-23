/**
 * This program test whether the state of Mersenne Twister
 * given in the folder mersenne_twister_states are correctly
 * separated by a fixed amount of random numbers in the sequence
 * of Mersenne Twister
 */

#include <stdio.h>
#include "../libmt/mt19937ar.h"
#define N 624

#define GAP 1000000000000 // amount of random number between 2 states
#define TMP_FILENAME "tmp"

/* compares the states stored in 2 files and returns a boolean */
int areStateFileEqual(char* fileA, char* fileB);

int main(int argc, char const *argv[])
{
    char* filenames[] = { "../mersenne_twister_states/mts000M000000", // first test
                          "../mersenne_twister_states/mts000M000001",
                          "../mersenne_twister_states/mts000M000202", // second test
                          "../mersenne_twister_states/mts000M000203" };
    int i;
    unsigned long j;

    // no initialisation because with restore an existing status

    for(i = 0; i < 4; i+=2)
    {
        restoreStatus(filenames[i]);

        // generate the amount of random numbers required to
        // reach the next state
        for(j = 0 ; j<GAP ; j++)
            genrand_int32();

        // compare the states
        saveStatus(TMP_FILENAME);
        if( !areStateFileEqual(TMP_FILENAME, filenames[i+1]) )
        {
            return 1; // test failed
        }
    }

    return 0; // test succeed
}

int areStateFileEqual(char* fileA, char* fileB)
{
    int mtiA, mtiB, i;
    unsigned long mtA, mtB;

    FILE * fA  = fopen(fileA, "r");
    FILE * fB  = fopen(fileB, "r");
    if (fA == NULL || fB == NULL)
        return 0; // false, they are not equal if one of them can't be openned

    fscanf(fA, "%d", &mtiA);
    fscanf(fB, "%d", &mtiB);
    if(mtiA != mtiB)
    {
        fclose(fA);
        fclose(fB);
        return 0; // false
    }

    
    for (i = 0; i < N; i++)
    {
        fscanf(fA, "%ld", &mtA);
        fscanf(fB, "%ld", &mtB);
        if (mtA != mtB)
        {
            fclose(fA);
            fclose(fB);
            return 0; // false
        }
    }
    
    fclose(fA);
    fclose(fB);

    return 1; // true
}