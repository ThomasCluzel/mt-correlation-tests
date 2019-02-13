#include <iostream>
#include "allexperiments.hpp"

const int MAX_NUMBER_OF_STATUS = 128;

using namespace std;

void runExperiments(const char* filesWithStatus[MAX_NUMBER_OF_STATUS], int flagTestsEnable, double alpha)
{
    for(int numberOfStatus=2; numberOfStatus<=MAX_NUMBER_OF_STATUS; numberOfStatus*=2) // for each (n, p) possible
    {
        cout << "\n\nExperiments with " << numberOfStatus << " status\n\n";
        runCorrelationTests(filesWithStatus, numberOfStatus, flagTestsEnable, alpha); // run tests
    }
}

void runControlExperiments()
{
    // initialize the status
    const char* initialStatus = "mersenne_twister_states/mts000M000000";
    const char* statusFiles[MAX_NUMBER_OF_STATUS];
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
        statusFiles[i] = initialStatus;

    cout << "\nBeginning of the control experiments (all status are the same)\n";

    // run the control experiments
    runExperiments(statusFiles, FLAG_TEST_CORR | FLAG_TEST_MULT);
}

void runExperimentsStatusSeparatedByOneTrillionNumbers()
{
    // each status are separated by one trillion number
    const char* baseName = "mersenne_twister_states/mts000M000";
    char* statusFiles[MAX_NUMBER_OF_STATUS]; // array of pointers, not array of arrays
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
    {
        statusFiles[i] = (char*)malloc(38 * sizeof(char));
        sprintf(statusFiles[i], "%s%03d", baseName, i);
    }

    cout << "\nBeginning of the experiments with status separated by one trillion numbers\n";

    // run the all the experiments
    runExperiments((const char**)statusFiles);

    // free memory
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
        free(statusFiles[i]);
}