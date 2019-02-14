#include <iostream>
#include "allexperiments.hpp"

const int MAX_NUMBER_OF_STATUS = 128;

using std::cout;

void runExperiments(const char* filesWithStatus[MAX_NUMBER_OF_STATUS], int flagTestsEnable, double alpha)
{
    for(int numberOfStatus=2; numberOfStatus<=MAX_NUMBER_OF_STATUS; numberOfStatus*=2) // for each (n, p) possible
    {
        cout << "\n\n\nExperiments with " << numberOfStatus << " status\n\n";
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

    // run the control experiments
    runExperiments(statusFiles, FLAG_TEST_CORR | FLAG_TEST_MULT);
}

void runAllExperiments(int step)
{
    // each status are separated by one trillion number
    const char* baseName = "mersenne_twister_states/mts000M000";
    char* statusFiles[MAX_NUMBER_OF_STATUS]; // array of pointers, not array of arrays
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i+=step)
    {
        statusFiles[i] = (char*)malloc(38 * sizeof(char));
        sprintf(statusFiles[i], "%s%03d", baseName, i);
    }

    // run the all the experiments
    runExperiments((const char**)statusFiles);

    // free memory
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i+=step)
        free(statusFiles[i]);
}